#!/usr/bin/env luajit
local cmdline = require 'ext.cmdline'(...)
local assertindex = require 'ext.assert'.index
local gl = require 'gl.setup'(cmdline.gl or 'OpenGL')
local ffi = require 'ffi'
local vec3f = require 'vec-ffi.vec3f'
local vec4f = require 'vec-ffi.vec4f'
local quatf = require 'vec-ffi.quatf'
local matrix_ffi = require 'matrix.ffi'
local template = require 'template'
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local path = require 'ext.path'
local glreport = require 'gl.report'
local ig = require 'imgui'
local GLTex2D = require 'gl.tex2d'
local GLFBO = require 'gl.fbo'
local GLArrayBuffer = require 'gl.arraybuffer'
local GLProgram = require 'gl.program'
local GLGeometry = require 'gl.geometry'
local GLSceneObject = require 'gl.sceneobject'
local clnumber = require 'cl.obj.number'
local StatSet = require 'stat.set'


-- Sets WGS-84 parameters
local wgs84 = {}
wgs84.a = 6378.137 --semi-major axis of the ellipsoid in
wgs84.b = 6356.7523142 --semi-minor axis of the ellipsoid in
wgs84.fla = 1 / 298.257223563 -- flattening
wgs84.eps = math.sqrt(1 - (wgs84.b * wgs84.b) / (wgs84.a * wgs84.a)) --first eccentricity
wgs84.epssq = wgs84.eps * wgs84.eps --first eccentricity squared
wgs84.re = 6371.2 -- Earth's radius


local HeightAboveEllipsoid = 0		-- compute at z=0 for now
-- TODO add support for dg/dh

-- load wmm
local lines = string.split(string.trim(path'wmm.cof':read()), '\n')
local header = string.split(string.trim(lines:remove(1)), '%s+')
assert(#header == 3)
local wmm = {}
wmm.epoch = tonumber(header[1])
print('model epoch', wmm.epoch)
print('model name', header[2])
print('date of release', header[3])
assert(lines:remove() == '999999999999999999999999999999999999999999999999')
assert(lines:remove() == '999999999999999999999999999999999999999999999999')

for _,line in ipairs(lines) do
	local parts = string.split(string.trim(line), '%s+')
	assert(#parts == 6)
	parts = parts:mapi(function(x)
		return assert(tonumber(x))
	end)
	local n = parts[1]
	local m = parts[2]
	wmm[n] = wmm[n] or {}
	wmm[n][m] = wmm[n][m] or {}
	wmm[n][m] = {
		g = parts[3],
		h = parts[4],
		dg_dt = parts[5],
		dh_dt = parts[6],
	}
	-- n runs 1-12
	-- m runs 0-n
end

local nMax = #wmm

local function calcB(phi, lambda, height)

	-- begin MAG_GeodeticToSpherical
	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)

	-- convert from geodetic WGS-84 to spherical coordiantes
	local rc = wgs84.a / math.sqrt(1 - wgs84.epssq * sinPhi * sinPhi)
	local xp = (rc + height) * cosPhi
	local zp = (rc * (1 - wgs84.epssq) + height) * sinPhi

	-- spherical results:
	local r = math.sqrt(xp * xp + zp * zp)
	local sinPhiSph = zp / r	-- geocentric latitude sin & cos
	local cosPhiSph = math.sqrt(1 - sinPhiSph * sinPhiSph)
	-- longitude is the same
	-- end MAG_GeodeticToSpherical

	-- begin MAG_Geomag
	-- begin MAG_ComputeSphericalHarmonicVariables

	local cosLambda = math.cos(lambda)
	local sinLambda = math.sin(lambda)

	local RelativeRadiusPower = {}	-- 0-nmake
	RelativeRadiusPower[0] = (wgs84.re / r) * (wgs84.re / r)
	for n=1,nMax do
		RelativeRadiusPower[n] = RelativeRadiusPower[n-1] * (wgs84.re / r)
	end

	local cosLambdaToTheM = {[0] = 1, cosLambda}
	local sinLambdaToTheM = {[0] = 0, sinLambda}

	-- looks like exp(i lambda)
	for m=2,nMax do
		cosLambdaToTheM[m] = cosLambdaToTheM[m-1] * cosLambda - sinLambdaToTheM[m-1] * sinLambda
		sinLambdaToTheM[m] = cosLambdaToTheM[m-1] * sinLambda + sinLambdaToTheM[m-1] * cosLambda
	end

	-- end MAG_ComputeSphericalHarmonicVariables
	-- begin MAG_AssociatedLegendreFunction

	-- begin MAG_PcupLow

	local x = sinPhiSph
	local Pcup = {[0] = 1}
	local dPcup = {[0] = 0}

	-- sin (geocentric latitude) - sinPhiSph
	local z = math.sqrt((1 - x) * (1 + x));

	local schmidtQuasiNorm = {}

	--	 First,	Compute the Gauss-normalized associated Legendre  functions
	for n=1,nMax do
		for m=0,n do
			local index = n * (n + 1) / 2 + m
			if n == m then
				local index1 = (n - 1) * n / 2 + m - 1
				Pcup[index] = z * Pcup[index1]
				dPcup[index] = z * dPcup[index1] + x * Pcup[index1]
			elseif n == 1 and m == 0 then
				local index1 = (n - 1) * n / 2 + m
				Pcup[index] = x * Pcup[index1]
				dPcup[index] = x * dPcup[index1] - z * Pcup[index1]
			elseif n > 1 and n ~= m then
				local index1 = (n - 2) * (n - 1) / 2 + m
				local index2 = (n - 1) * n / 2 + m
				if m > n - 2 then
					Pcup[index] = x * Pcup[index2]
					dPcup[index] = x * dPcup[index2] - z * Pcup[index2]
				else
					local k = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3))
					Pcup[index] = x * Pcup[index2] - k * Pcup[index1]
					dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1]
				end
			end
		end
	end
	-- Compute the ration between the the Schmidt quasi-normalized associated Legendre
	-- functions and the Gauss-normalized version. */

	schmidtQuasiNorm[0] = 1;
	for n=1,nMax do
		local index = (n * (n + 1) / 2);
		local index1 = (n - 1) * n / 2;
		-- for m = 0
		schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (2 * n - 1) / n

		for m=1,n do
			local index = (n * (n + 1) / 2 + m)
			local index1 = (n * (n + 1) / 2 + m - 1)
			schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * math.sqrt( ((n - m + 1) * (m == 1 and 2 or 1)) / (n + m));
		end
	end

	-- Converts the  Gauss-normalized associated Legendre
	-- functions to the Schmidt quasi-normalized version using pre-computed
	-- relation stored in the variable schmidtQuasiNorm */

	for n=1,nMax do
		for m=0,n do
			local index = (n * (n + 1) / 2 + m)
			Pcup[index] = Pcup[index] * schmidtQuasiNorm[index]
			dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index]
			-- The sign is changed since the new WMM routines use derivative with respect to latitude
			-- insted of co-latitude */
		end
	end

	-- end MAG_PcupLow
	-- end MAG_AssociatedLegendreFunction
	-- begin MAG_Summation

	local Bz = 0
	local By = 0
	local Bx = 0
	for n=1,nMax do
		local wmm_n = wmm[n]
		for m=0,n do
			local index = (n * (n + 1) / 2 + m)
			local wmm_n_m = wmm_n[m]

			--		    nMax  	(n+2) 	  n     m            m           m
			--		Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
			--						n=1      	      m=0   n            n           n  */
			-- Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
			Bz = Bz -
				RelativeRadiusPower[n]
				* (
					wmm_n_m.g * cosLambdaToTheM[m]
					+ wmm_n_m.h * sinLambdaToTheM[m]
				)
				* (n + 1) * Pcup[index]

			--		  1 nMax  (n+2)    n     m            m           m
			--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1             m=0   n            n           n  */
			-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
			By = By + (
				RelativeRadiusPower[n]
				* (
					wmm_n_m.g * sinLambdaToTheM[m]
					- wmm_n_m.h * cosLambdaToTheM[m]
				)
				* m * Pcup[index]
			)
			--		   nMax  (n+2) n     m            m           m
			--		Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1         m=0   n            n           n  */
			-- Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */
			Bx = Bx -
				RelativeRadiusPower[n] *
				(
					wmm_n_m.g * cosLambdaToTheM[m]
					+ wmm_n_m.h * sinLambdaToTheM[m]
				)
				* dPcup[index]
		end
	end

	if cosPhiSph < -1e-10 or cosPhiSph > 1e-10 then
		By = By / cosPhiSph
	else
		-- Special calculation for component - By - at Geographic poles.
		-- If the user wants to avoid using this function,  please make sure that
		-- the latitude is not exactly +/-90. An option is to make use the function
		-- MAG_CheckGeographicPoles.
		-- begin MAG_SummationSpecial

		local PcupS = {[0] = 1}
		local schmidtQuasiNorm1 = 1

		By = 0

		for n=1,nMax do
			--Compute the ration between the Gauss-normalized associated Legendre
			-- functions and the Schmidt quasi-normalized version. This is equivalent to
			-- sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)!  */
			local m = 1
			local wmm_n_m = wmm[n][m]

			local index = (n * (n + 1) / 2 + m)
			local schmidtQuasiNorm2 = schmidtQuasiNorm1 *  (2 * n - 1) /  n
			local schmidtQuasiNorm3 = schmidtQuasiNorm2 * math.sqrt( (n * 2) /  (n + 1))
			local schmidtQuasiNorm1 = schmidtQuasiNorm2
			if n == 1 then
				PcupS[n] = PcupS[n-1]
			else
				local k =  (((n - 1) * (n - 1)) - 1) /  ((2 * n - 1) * (2 * n - 3))
				PcupS[n] = sinPhiSph * PcupS[n - 1] - k * PcupS[n - 2]
			end

			--		  1 nMax  (n+2)    n     m            m           m
			--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1             m=0   n            n           n  */
			-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
			By = By +
				RelativeRadiusPower[n] *
				(
					wmm_n_m.g * sinLambdaToTheM[1] -
					wmm_n_m.h * cosLambdaToTheM[1])
					* PcupS[n] * schmidtQuasiNorm3
		end

		-- end MAG_SummationSpecial
	end

	-- end MAG_Summation
	-- end MAG_Geomag

	return Bx, By, Bz
end

-- expects xyz in cartesian units of wgs84's 'a' parameter
-- TODO just use charts.WGS84:chartInv(x,y,z) ?
local function cartesianToLatLonWGS84(x, y, z)
	x = x * wgs84.a
	y = y * wgs84.a
	z = z * wgs84.a

	local modified_b = z < 0 and -wgs84.b or wgs84.b

	local r = math.sqrt(x*x + y*y)

	local e = ( modified_b*z - (wgs84.a*wgs84.a - modified_b*modified_b) ) / ( wgs84.a*r )
	local f = ( modified_b*z + (wgs84.a*wgs84.a - modified_b*modified_b) ) / ( wgs84.a*r )
	local p = (4 / 3) * (e*f + 1)
	local q = 2 * (e*e - f*f)
	local d = p*p*p + q*q

	local v
	if  d >= 0 then
		v = math.pow( (math.sqrt( d ) - q), (1 / 3) )
			- math.pow( (math.sqrt( d ) + q), (1 / 3) )
	else
		v= 2 * math.sqrt( -p )
			* math.cos( math.acos( q/(p * math.sqrt( -p )) ) / 3 )
	end

	if v*v < math.abs(p)  then
		v = -(v*v*v + 2*q) / (3*p)
	end

	local g = (math.sqrt( e*e + v ) + e) / 2
	local t = math.sqrt( g*g  + (f - v*g)/(2*g - e) ) - g

	local rlat = math.atan( (wgs84.a*(1 - t*t)) / (2*modified_b*t) )
	local phi = rlat

	local height = (r - wgs84.a*t) * math.cos(rlat) + (z - modified_b) * math.sin(rlat)
	local zlong = math.atan2(y, x)
	if  zlong < 0 then
		zlong = zlong + 2*math.pi
	end
	local lambda = zlong
	while lambda > math.pi do
		lambda= lambda - 2 * math.pi
	end
	return phi, lambda, height
end

--[[
Bx = d/dphi points north
By = d/dlambda points east
Bz = -d/dheight points inwards
--]]

local BStat = StatSet(
	'mag', 'x', 'y', 'z', 'mag2d',
	'div', 'div2d', 'curlZ', 'curlMag'
)
for i,stat in ipairs(BStat) do stat.onlyFinite = true end

local londim -- dimension in 2D x dir / spherical phi / globe lambda dir
local latdim -- dimension in 2D y dir / spherical theta / globe phi dir

local BTex
local B2Tex
local earthtex

local App = require 'imguiapp.withorbit'()
App.viewUseBuiltinMatrixMath = true
App.title = 'EM field'


local guivars = {
	geomIndex = 1,
	overlayIndex = 1,
	gradientIndex = 0,

	drawAlpha = 1,
	doDrawVectorField = true,

	fieldDT = 0,


	-- set to 0 to flatten vector field against surface
	fieldZScale = 1,

	arrowScale = 5,

	gradScale = 1,
}



local Geom = class()

function Geom:init(args)
	for k,v in pairs(args) do
		self[k] = v
	end
end

local allCharts = require 'geographic-charts'	-- TODO switch to this

-- which charts we want to allow ...
local chartNames = table{
	'Equirectangular',
	'Azimuthal equidistant',
	'Mollweide',
	'WGS84',
}

local chartIndexForName = chartNames:mapi(function(v,i) return i,v end)

local charts = table()
for i,name in ipairs(chartNames) do
	local chart = assertindex(allCharts, name)
	charts[i] = chart
	charts[name] = chart
end

charts['Azimuthal equidistant'].R = .5 * math.pi
--allCharts.WGS84_a = wgs84.a
for _,name in ipairs(chartNames) do
	local c = charts[name]
	if c.build then c:build() end
end

local chartCNames = charts:mapi(function(chart) return chart:getCName() end)

-- TODO this but only for the charts we're using ...
--require 'geographic-charts.buildall'
local chartCode = require 'geographic-charts.code'(charts)
--local chartCode = modules:getCodeAndHeader(table.mapi(charts, function(c) return c:getSymbol() end):unpack())

for _,c in ipairs(charts) do
	local oldChartFunc = c.chart
	function c:chart(phi, lambda, height)
		local x, y, z = oldChartFunc(self, math.deg(phi), math.deg(lambda), height)
		x = x / allCharts.WGS84_a
		y = y / allCharts.WGS84_a
		z = z / allCharts.WGS84_a
		return x, y, z
	end

	local oldBasisFunc = c.basis
	function c:basis(phi, lambda, height)
		local ex, ey, ez = oldBasisFunc(self, math.deg(phi), math.deg(lambda), height)
		return ex, ey, ez
	end
end

local gradients = {
	{
		name = 'rainbow',
		gen = function()
			local tex = require 'gl.hsvtex2d'(256)
				:setWrap{
					s = gl.GL_REPEAT,
					t = gl.GL_REPEAT,
				}
				:unbind()
			return tex
		end,
	},
	{
		name = 'stripes',
		gen = function()
			local n = 256
			-- why is this intermittantly misaligned?
			local image = require 'image'(
				n, 1, 4, 'unsigned char', function(u)
					local s = (u+.5)/n
					-- TODO this but in shader so we can dynamically change it
					local alpha = bit.band(u, 15) == 0 and 255 or 0
					return 255*(1-s),0,255*s,alpha
				end
			)
			return GLTex2D{
				--image = image,
				data = image.buffer,
				width = image.width,
				height = 1,

				internalFormat = gl.GL_RGBA,
				format = gl.GL_RGBA,
				type = gl.GL_UNSIGNED_BYTE,
				minFilter = gl.GL_LINEAR,
				--magFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
				magFilter = gl.GL_LINEAR,
				wrap = {
					s = gl.GL_REPEAT,
					t = gl.GL_REPEAT,
				},
				generateMipmap = true,
			}:unbind()
		end,
	},
}

local overlays = {
	{
		name = 'Earth',
		code = [[
	hsvBlend = 0.;
]],
	},
	{
		name = '|B|',
		code = [[
	s = (length(B) - BMagMin) / (BMagMax - BMagMin);
]],
	},
	{
		name = '|B| 2D',
		code = [[
	s = (length(B.xy) - BMag2DMin) / (BMag2DMax - BMag2DMin);
]],
	},
	{
		name = 'arg B',
		code = [[
	s = atan(B.y, B.x) / (2. * M_PI) + .5;
]],
	},
	{
		name = 'Bx (north)',
		code = [[
	s = (B.x - BxMin) / (BxMax - BxMin);
]],
	},
	{
		name = 'By (east)',
		code = [[
	s = (B.y - ByMin) / (ByMax - ByMin);
]],
	},
	{
		name = 'Bz (inward)',
		code = [[
	s = (B.z - BzMin) / (BzMax - BzMin);
]],
	},
	{
		name = 'div B',
		code = [[
	s = (B2.x - BDivMin) / (BDivMax - BDivMin);
]],
	},
	{
		name = 'div2D B',
		code = [[
	s = (B2.y - BDiv2DMin) / (BDiv2DMax - BDiv2DMin);
]],
	},
	{
		name = 'curl B z',
		code = [[
	s = (B2.z - BCurlZMin) / (BCurlZMax - BCurlZMin);
]],
	},
	{
		name = 'curl B Mag',
		code = [[
	s = (B2.w - BCurlMagMin) / (BCurlMagMax - BCurlMagMin);
]],
	},
}

for chartIndex,c in ipairs(charts) do
	function c:draw(app, shader, gradtex)
		local height = 0
		local jres = 360
		local ires = 180
		if not self.sceneobj then
			local vertexes = table()
			for i=0,ires do
				local u = i/ires
				local phi = math.rad((u * 2 - 1) * 90)
				for j=0,jres do
					local v = j/jres
					local lambda = math.rad((v * 2 - 1) * 180)
					vertexes:append{v, u}	-- lon, lat = u, v in texcoord space
				end
			end
			self.vertexBuf = GLArrayBuffer{
				data = vertexes,
				dim = 2,
			}:unbind()

			local geometries = table()
			for ibase=0,ires-1 do
				local indexes = table()
				for j=0,jres do
					for iofs=1,0,-1 do
						local i = ibase + iofs
						indexes:insert(j + (jres + 1) * i)
					end
				end
				geometries:insert(GLGeometry{
					mode = gl.GL_TRIANGLE_STRIP,
					indexes = {
						data = indexes,
					},
					vertexes = self.vertexBuf,
				})
			end

			self.sceneobj = GLSceneObject{
				program = shader,
				vertexes = self.vertexBuf,
				geometries = geometries,
				texs = {earthtex, gradtex, BTex, B2Tex},
			}
		end

		self.sceneobj.program = shader
		self.sceneobj.texs[2] = gradtex
		self.sceneobj.uniforms.mvProjMat = app.view.mvProjMat.ptr
		self.sceneobj.uniforms['weight_Equirectangular'] = chartIndex == 1 and 1 or 0
		self.sceneobj.uniforms['weight_Azimuthal_equidistant'] = chartIndex == 2 and 1 or 0
		self.sceneobj.uniforms['weight_Mollweide'] = chartIndex == 3 and 1 or 0
		self.sceneobj.uniforms['weight_WGS84'] = chartIndex == 4 and 1 or 0
		self.sceneobj.uniforms.chartIs3D = self.is3D or false
		self.sceneobj:draw()
	end
end


function App:initGL(...)
	if App.super.initGL then
		App.super.initGL(self, ...)
	end

	earthtex = GLTex2D{
		filename = 'earth.png',
		-- hmm, default magfilter on my amd ubuntu seems to expect mipmaps, but doesn't generate mipmaps...?
		-- so generate them too
		generateMipmap = true,
	}

	for _,grad in ipairs(gradients) do
		grad.tex = grad:gen()
	end

	local calc_b_shader = template(assert(path'calc_b.shader':read()), {
		wgs84 = wgs84,
		wmm = wmm,
	})

	self.quadGeom = GLGeometry{
		mode = gl.GL_TRIANGLE_STRIP,
		vertexes = {
			data = {
				0, 0,
				1, 0,
				0, 1,
				1, 1,
			},
			dim = 2,
		},
	}

	self.unitProjMat = matrix_ffi({4,4}, 'float'):zeros():setOrtho(0, 1, 0, 1, -1, 1)
	--self.unitProjMat = matrix_ffi({4,4}, 'float'):zeros():setOrtho(-1, 1, -1, 1, 1, -1)	-- identity matrix

	do
		print'generating B field...'

		-- can be arbitrary
		-- but the WMM model is for 15 mins, so [360,180] x4
		londim = 1440
		latdim = 720

		-- while we're here, generate the basis tex
		-- i've only got a script implementation for it for now
		-- but i could generate a GLSL one from the expressions ...
--[[
require 'vec-ffi.quatf'
		local basisTexData = ffi.new('quatf_t[?]', londim * latdim)
print"calculating basis..."
		for _,chart in ipairs(charts) do
			local index = 0
			for j=0,latdim-1 do
				local v = (j + .5) / latdim
				local lat = (v * 2 - 1) * 90
				local phi = math.rad(lat)

				for i=0,londim-1 do
					local u = (i + .5) / londim
					local lon = (u * 2 - 1) * 180
					local lambda = math.rad(lon)

					local ex, ey, ez = chart:basis(phi, lambda, 0)
					ez = -ez	-- left-hand to right-hand
--DEBUG:require 'ext.assert'.ge(ex:cross(ey):dot(ez), .9)
					basisTexData[index]:fromMatrix{ex, ey, ez}:normalize()
					index = index + 1
				end
			end

			chart.basisTex = GLTex2D{
				internalFormat = gl.GL_RGBA32F,
				width = londim,
				height = latdim,
				format = gl.GL_RGBA,
				type = gl.GL_FLOAT,
				minFilter = gl.GL_LINEAR,
				magFilter = gl.GL_NEAREST,
				wrap = {
					s = gl.GL_REPEAT,
					t = gl.GL_REPEAT,
				},
				data = basisTexData,
				generateMipmap = true,
			}:unbind()
			chart.basisTex.data = nil	-- don't keep track of it since we're reusing this buffer
		end
print"...done calculating basis"
glreport'here'
--]]

		local fbo = GLFBO()
			:unbind()
glreport'here'

		local quadGeomVertexCode = [[
layout(location=0) in vec2 vertex;
out vec2 texcoordv;
uniform mat4 mvProjMat;
void main() {
	texcoordv = vertex;
	gl_Position = mvProjMat * vec4(vertex, 0., 1.);
}
]]

		local calcBShader = GLProgram{
			version = 'latest',
			precision = 'best',
			vertexCode = quadGeomVertexCode,
			fragmentCode = [[
uniform float dt;
]]..calc_b_shader..template[[

#define M_PI <?=('%.49f'):format(math.pi)?>

in vec2 texcoordv;
out vec4 fragColor;
void main() {
	float phi = (texcoordv.y - .5) * M_PI;			//[-pi/2, pi/2]
	float lambda = (texcoordv.x - .5) * 2. * M_PI;	//[-pi, pi]
	vec3 B = calcB(vec3(phi, lambda, 0.));
	fragColor = vec4(B, 1.);
}
]],
			uniforms = {
				dt = 0,
			},
		}:useNone()
glreport'here'

		BTex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,
			width = londim,
			height = latdim,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
			magFilter = gl.GL_NEAREST,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}:unbind()
glreport'here'

		local calcBSceneObj = GLSceneObject{
			program = calcBShader,
			geometry = self.quadGeom,
		}

		-- BData / B2Data is only used for stat computation
		local BData = ffi.new('vec4f_t[?]', londim * latdim)
		fbo:draw{
			viewport = {0, 0, londim, latdim},
			dest = BTex,
			callback = function()
				gl.glClear(gl.GL_COLOR_BUFFER_BIT)
				calcBSceneObj.uniforms.mvProjMat = self.unitProjMat.ptr
				calcBSceneObj:draw()
				gl.glReadPixels(0, 0, BTex.width, BTex.height, gl.GL_RGBA, gl.GL_FLOAT, BData)
			end,
		}
glreport'here'
		BTex:bind()
			:generateMipmap()
			:unbind()
glreport'here'


-- [=[ hmm, better way than copy paste?

		print'generating div B and curl B...'

		local calcB2Shader = GLProgram{
			version = 'latest',
			precision = 'best',
			vertexCode = quadGeomVertexCode,
			fragmentCode = [[
uniform float dt;
]]..calc_b_shader..template([[

#define M_PI <?=('%.49f'):format(math.pi)?>

#define latdim	<?=clnumber(latdim)?>
#define londim	<?=clnumber(londim)?>

vec3 dphi = vec3(M_PI / londim, 0., 0.);
vec3 dlambda = vec3(0., 2. * M_PI / latdim, 0.);
vec3 dheight = vec3(0., 0., 1.);

in vec2 texcoordv;
out vec4 fragColor;

void main() {
	float phi = (texcoordv.y - .5) * M_PI;			//[-pi/2, pi/2]
	float lambda = (texcoordv.x - .5) * 2. * M_PI;	//[-pi, pi]

	vec3 plh = vec3(phi, lambda, 0.);

	vec3 dphi_B = (calcB(plh + dphi) - calcB(plh - dphi)) / dphi.x / (wgs84_a * cos(plh.x));
	vec3 dlambda_B = (calcB(plh + dlambda) - calcB(plh - dlambda)) / dlambda.y / wgs84_a;
	vec3 dheight_B = (calcB(plh + dheight) - calcB(plh - dheight)) / dheight.z;

	float div2D_B = dphi_B.x + dlambda_B.y;
	float div_B = div2D_B + dheight_B.z;

	vec3 curl_B = vec3(
		dlambda_B.z - dheight_B.y,
		dheight_B.x - dphi_B.z,
		dphi_B.y - dlambda_B.x
	);

	fragColor = vec4(
		div_B,
		div2D_B,
		curl_B.z,
		length(curl_B)
	);
}
]],
				{
					latdim = latdim,
					londim = londim,
					clnumber = clnumber,
				}
			),
			uniforms = {
				dt = 0,
			},
		}:useNone()
glreport'here'

		B2Tex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,
			width = londim,
			height = latdim,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
			magFilter = gl.GL_NEAREST,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}:unbind()
glreport'here'

		local calcB2SceneObj = GLSceneObject{
			program = calcB2Shader,
			geometry = self.quadGeom,
		}

		-- only used for stat calc
		local B2Data = ffi.new('vec4f_t[?]', londim * latdim)
		fbo:draw{
			viewport = {0, 0, londim, latdim},
			dest = B2Tex,
			callback = function()
				gl.glClear(gl.GL_COLOR_BUFFER_BIT)
				calcB2SceneObj.uniforms.mvProjMat = self.unitProjMat.ptr
				calcB2SceneObj:draw()
				gl.glReadPixels(0, 0, B2Tex.width, B2Tex.height, gl.GL_RGBA, gl.GL_FLOAT, B2Data)
			end,
		}
glreport'here'
		B2Tex:bind()
			:generateMipmap()
			:unbind()

		local statgens = table{
			function(B, B2) return math.sqrt(B.x*B.x + B.y*B.y + B.z*B.z) end,	-- not :length() since it's vec4...
			function(B, B2) return B.x end,
			function(B, B2) return B.y end,
			function(B, B2) return B.z end,
			function(B, B2) return math.sqrt(B.x*B.x + B.y*B.y) end,
			function(B, B2) return B2.x end,
			function(B, B2) return B2.y end,
			function(B, B2) return B2.z end,
			function(B, B2) return B2.w end,
		}

		for j=0,latdim-1 do
			for i=0,londim-1 do
				local e = i + londim * j
				local B = BData[e]
				local B2 = B2Data[e]
				local Bmag = B:length()
				BStat:accum(
					statgens:mapi(function(f)
						return f(B, B2)
					end):unpack()
				)
			end
		end

		print('BStat')
		print(BStat)

--[==[	-- plotting the bins
		local Bin = require 'stat.bin'

		local binCount = 100
		local bins = table.mapi(BStat, function(stat,k)
			return Bin(
				math.max(stat.min, stat.avg - 3*stat.stddev),
				math.min(stat.max, stat.avg + 3*stat.stddev),
				binCount
			)
		end)
		for j=0,latdim-1 do
			for i=0,londim-1 do
				local e = i + londim * j
				local B = BData[e]
				local B2 = B2Data[e]
				local Bmag = B:length()
				for i=1,#bins do
					local stat = BStat[i]
					local v = statgens[i](B, B2)
					if v >= stat.avg - 3*stat.stddev and v <= stat.avg + 3*stat.stddev then
						bins[i]:accum(v)
					end
				end
			end
		end
		for i,bin in ipairs(bins) do
			path'tmp.txt':write(bin:getTextData())
			require 'gnuplot'{
				title = BStat[i].name,
				terminal = 'png size 1024,768',
				output = BStat[i].name..'.png',
				{using='1:2', datafile = 'tmp.txt'},
			}
		end
--]==]

		-- clamp in the min/max to 3 stddev
		for i=1,#BStat do
			local stat = BStat[i]
			stat.min = math.max(stat.min, stat.avg - 3*stat.stddev)
			stat.max = math.min(stat.max, stat.avg + 3*stat.stddev)
		end
--]=]
	end

	glreport'here'

	local overlayVertexCode = chartCode..template([[

<? for _,name in ipairs(chartCNames) do
?>uniform float weight_<?=name?>;
<? end
?>
uniform bool chartIs3D;

layout(location=0) in vec2 vertex;
out vec2 texcoordv;
uniform mat4 mvProjMat;
void main() {
	texcoordv = vertex;	//(lon, lat) in [0,1]

	vec3 coords = vec3(
		(vertex.y - .5) * 180.,	// lat in deg, [-90, 90]
		(vertex.x - .5) * 360.,	// lon in deg, [-180, 180]
		0.);						// height in meters

	// expect vertex xyz to be lat lon height
	// lat and lon is in degrees
	// height is in meters
	vec3 pos = 0.
<? for _,name in ipairs(chartCNames) do
?>		+ weight_<?=name?> * chart_<?=name?>(coords)
<? end
?>	;
	//from meters to normalized coordinates
	pos /= WGS84_a;

	if (chartIs3D) {
		pos = vec3(pos.z, pos.x, pos.y);
	}

	gl_Position = mvProjMat * vec4(pos, 1.);
}
]],
	{
		chartCNames = chartCNames,
	})

	for i,overlay in ipairs(overlays) do
		overlay.shader = GLProgram{
			version = 'latest',
			precision = 'best',
			vertexCode = overlayVertexCode,
			fragmentCode = [[
uniform float dt;
]]..calc_b_shader..template([[

// TODO if you use tex then at the very pole there is a numeric artifact ...
#define CALC_B_ON_GPU

#define BMagMin <?=clnumber(BStat.mag.min)?>
#define BMagMax <?=clnumber(BStat.mag.max)?>
#define BMag2DMin <?=clnumber(BStat.mag2d.min)?>
#define BMag2DMax <?=clnumber(BStat.mag2d.max)?>
#define BxMin <?=clnumber(BStat.x.min)?>
#define BxMax <?=clnumber(BStat.x.max)?>
#define ByMin <?=clnumber(BStat.y.min)?>
#define ByMax <?=clnumber(BStat.y.max)?>
#define BzMin <?=clnumber(BStat.z.min)?>
#define BzMax <?=clnumber(BStat.z.max)?>
#define BDivMin <?=clnumber(BStat.div.min)?>
#define BDivMax <?=clnumber(BStat.div.max)?>
#define BDiv2DMin <?=clnumber(BStat.div2d.min)?>
#define BDiv2DMax <?=clnumber(BStat.div2d.max)?>
#define BCurlZMin <?=clnumber(BStat.curlZ.min)?>
#define BCurlZMax <?=clnumber(BStat.curlZ.max)?>
#define BCurlMagMin <?=clnumber(BStat.curlMag.min)?>
#define BCurlMagMax <?=clnumber(BStat.curlMag.max)?>

#define M_PI <?=('%.49f'):format(math.pi)?>

in vec2 texcoordv;
out vec4 fragColor;

uniform sampler2D earthTex;
uniform sampler2D BTex;
uniform sampler2D B2Tex;
uniform sampler2D gradTex;
uniform float alpha;
uniform float gradScale;

void main() {
	float s = .5;
	float hsvBlend = .5;

	vec4 B2 = texture(B2Tex, texcoordv);
#ifndef CALC_B_ON_GPU
	vec4 B = texture(BTex, texcoordv);
#else
	vec3 plh = vec3(
		(texcoordv.y - .5) * M_PI,			//phi
		(texcoordv.x - .5) * 2. * M_PI,		//lambda
		0.
	);
	vec3 B = calcB(plh);
#endif

	<?=overlay.code?>

	fragColor = mix(
		texture(earthTex, vec2(texcoordv.x, 1. - texcoordv.y)),
		texture(gradTex, vec2(s * gradScale, .5)),
		hsvBlend);
	fragColor.a = alpha;
}
]],
				{
					overlay = overlay,
					BStat = BStat,
					clnumber = clnumber,
				}
			),
			uniforms = table({
				earthTex = 0,
				gradTex = 1,
				BTex = 2,
				B2Tex = 3,
				alpha = 1,
				dt = 0,
			}, chartCNames:mapi(function(name,j)
				return j==i and 1 or 0, 'weight_'..name
			end)):setmetatable(nil),
		}:useNone()
	end

	if self.view then
		self.view.ortho = true
		self.view.pos.z = 2
		self.view.orthoSize = 2
	end
	gl.glClearColor(0,0,0,0)


	self.pointSceneObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
in vec3 vertex;
uniform mat4 mvProjMat;
uniform float pointSize;
void main() {
	gl_Position = mvProjMat * vec4(vertex, 1.);
	gl_PointSize = pointSize;
}
]],
			fragmentCode = [[
uniform vec3 color;
out vec4 fragColor;
void main() {
	fragColor = vec4(color, 1.);
}
]],
		},
		vertexes = {
			dim = 3,
			count = 1,
			size = 3 * ffi.sizeof'float',
			type = gl.GL_FLOAT,
			data = ffi.new('vec3f_t[?]', 1),
		},
		geometry = {
			mode = gl.GL_POINTS,
		},
	}

	self.lineSceneObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
in vec3 vertex;
uniform mat4 mvProjMat;
void main() {
	gl_Position = mvProjMat * vec4(vertex, 1.);
}
]],
			fragmentCode = [[
uniform vec3 color;
out vec4 fragColor;
void main() {
	fragColor = vec4(color, 1.);
}
]],
		},
		vertexes = {
			dim = 3,
			count = 2,
			size = 2 * 3 * ffi.sizeof'float',
			type = gl.GL_FLOAT,
			data = ffi.new('vec3f_t[?]', 2),
		},
		geometry = {
			mode = gl.GL_LINES,
		},
	}

	self.vectorFieldShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = chartCode
..[[
uniform float dt;
]]..calc_b_shader
..template([[

<? for _,name in ipairs(chartCNames) do
?>uniform float weight_<?=name?>;
<? end
?>
uniform bool chartIs3D;

layout(location=0) in vec3 vertex;

in vec2 texcoord;

out vec3 colorv;

uniform mat4 mvProjMat;
uniform float arrowScale;
uniform float fieldZScale;

#if 0	//basis from texture
uniform sampler2D basisTex;
#else	//basis from attribs
in vec3 ex, ey, ez;
#endif

#if 0	// B coeffs using textures
uniform sampler2D BTex;
#endif

float lenSq(vec3 v) {
	return dot(v, v);
}

vec3 perpTo(vec3 v) {
	vec3 vx = vec3(0., -v.z, v.y);
	vec3 vy = vec3(v.z, 0., -v.x);
	vec3 vz = vec3(-v.y, v.x, 0.);
	vec3 lenSqs = vec3(lenSq(vx), lenSq(vy), lenSq(vz));
	if (lenSqs.x > lenSqs.y) {		// x > y
		if (lenSqs.x > lenSqs.z) {	// x > y, x > z
			return vx;
		} else {					// z > x > y
			return vz;
		}
	} else {						// y > x
		if (lenSqs.y > lenSqs.z) {	// y > x, y > z
			return vy;
		} else {					// z > y > x
			return vz;
		}
	}
}

vec3 xAxis(vec4 q) {
	return vec3(
		1. - 2. * (q.y * q.y + q.z * q.z),
		2. * (q.x * q.y + q.z * q.w),
		2. * (q.x * q.z - q.w * q.y));
}

vec3 yAxis(vec4 q) {
	return vec3(
		2. * (q.x * q.y - q.w * q.z),
		1. - 2. * (q.x * q.x + q.z * q.z),
		2. * (q.y * q.z + q.w * q.x));
}

vec3 zAxis(vec4 q) {
	return vec3(
		2. * (q.x * q.z + q.w * q.y),
		2. * (q.y * q.z - q.w * q.x),
		1. - 2. * (q.x * q.x + q.y * q.y));
}

void main() {
	vec3 coords = vec3(
		(texcoord.y - .5) * 180.,	// lat in deg
		(texcoord.x - .5) * 360., 	// lon in deg
		0.);						// height in m

	// expect vertex xyz to be lat lon height
	// lat and lon is in degrees
	// height is in meters
	vec3 pos = 0.
<? for _,name in ipairs(chartCNames) do
?>		+ weight_<?=name?> * chart_<?=name?>(coords)
<? end
?>	;
	//from meters to normalized coordinates
	pos /= WGS84_a;

	// TODO why correct for pos but not basis ?
	if (chartIs3D) {
		pos = vec3(pos.z, pos.x, pos.y);
	}

#if 0	// B coeffs using textures
	vec3 BCoeffs = texture(BTex, texcoord).xyz;
#elif 1	// on GPU
	vec3 BCoeffs = calcB(vec3(
		(texcoord.y - .5) * M_PI,			//phi
		(texcoord.x - .5) * 2. * M_PI,		//lambda
		0.
	));
#endif

	// BCoeffs.z points inwards
	// I just changed ez to consistently point outwards (maybe I shouldn't have)
	BCoeffs.z = -BCoeffs.z;

#if 0	//basis from texture ... seems very inaccurate ...
	vec4 basisQuat = normalize(texture(basisTex, texcoord));
	vec3 ex = xAxis(basisQuat);
	vec3 ey = yAxis(basisQuat);
	vec3 ez = -zAxis(basisQuat);
#endif

	vec3 B = arrowScale * (
		  ex * BCoeffs.x
		+ ey * BCoeffs.y
		+ ez * (BCoeffs.z * fieldZScale)
	);
	//vec3 Bey = perpTo(B);			// cross with furthest axis
	//vec3 Bey = vec3(-B.y, B.x, 0.);	// cross with z axis
	vec3 Bey = arrowScale * (
		  ey * BCoeffs.x
		- ex * BCoeffs.y
		+ ez * (BCoeffs.z * fieldZScale)
	);

	vec3 Bez = ez;

	gl_Position = mvProjMat * vec4(
		pos + vertex.x * B + vertex.y * Bey + vertex.z * Bez * fieldZScale,
		1.);
}
]],
		{
			chartCNames = chartCNames,
		}),
		fragmentCode = [[
out vec4 fragColor;
void main() {
	fragColor = vec4(1., 1., 1., 1.);
}
]],
		uniforms = {
			basisTex = 0,
			BTex = 1,
		},
	}:useNone()

	gl.glEnable(gl.GL_CULL_FACE)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
end

local function degmintofrac(deg, min, sec)
	return deg + (1/60) * (min + (1/60) * sec)
end

local function drawReading(info)
	local chart = info.chart

	local height = 1e-2
	local lat = degmintofrac(table.unpack(info.lat))
	local lon = degmintofrac(table.unpack(info.lon))


	local phi = math.rad(lat)
	local lambda = math.rad(lon)
	local headingRad = math.rad(info.heading)

	local x,y,z = chart:chart(phi, lambda, height)
	local ex, ey, ez = chart:basis(phi, lambda, height)
	local H = (ex * math.cos(headingRad) + ey * math.sin(headingRad)) * .2

	-- TODO chart should be transforms, and apply to all chart rendered
	self.pointSceneObj.uniforms.color = {1,1,0}
	self.pointSceneObj.uniforms.pointSize = 5
	self.pointSceneObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
	self.pointSceneObj.vertexes.buffer.data[0]:set(x,y,z)
	self.pointSceneObj.vertexes.buffer
		:bind()
		:updateData()
		:unbind()
	self.pointSceneObj:draw()

	self.lineSceneObj.uniforms.mvProjMat = self.view.mvProjMat.ptr

	self.lineSceneObj.uniforms.color = {0,0,1}
	self.lineSceneObj.vertexes.buffer.data[0]:set(x,y,z)
	self.lineSceneObj.vertexes.buffer.data[1]:set(x + ex.x * .2, y + ex.y * .2, z + ex.z * .2)
	self.lineSceneObj.vertexes.buffer
		:bind()
		:updateData()
		:unbind()
	self.lineSceneObj:draw()

	self.lineSceneObj.uniforms.color = {0,1,0}
	self.lineSceneObj.vertexes.buffer.data[0]:set(x,y,z)
	self.lineSceneObj.vertexes.buffer.data[1]:set(x + ey.x * .2, y + ey.y * .2, z + ey.z * .2)
	self.lineSceneObj.vertexes.buffer
		:bind()
		:updateData()
		:unbind()
	self.lineSceneObj:draw()

	self.lineSceneObj.uniforms.color = {0,0,1}
	self.lineSceneObj.vertexes.buffer.data[0]:set(x,y,z)
	self.lineSceneObj.vertexes.buffer.data[1]:set(x + H.x, y + H.y, z + H.z)
	self.lineSceneObj.vertexes.buffer
		:bind()
		:updateData()
		:unbind()
	self.lineSceneObj:draw()

	-- now use wgs84 here regardless of 'chart'
	--local pos = vec3f(latLonToCartesianWGS84(phi, lambda, height))
	local pos = vec3f(charts.WGS84:chart(math.rad(phi), math.rad(lambda), height)) / allCharts.WGS84_a
	local ex, ey, ez = charts.WGS84:basis(math.rad(phi), math.rad(lambda), .1)
	local H = (ex * math.cos(headingRad) + ey * math.sin(headingRad)) * .2

	local axis = pos:cross(H):normalize()

	local degToSpan = 90
	local numSteps = 90
	local dtheta = degToSpan / numSteps
	local q = quatf():fromAngleAxis(axis.x, axis.y, axis.z, dtheta)

	-- draw north vs magnetic north heading
	gl.glColor3f(1,1,1)
	gl.glBegin(gl.GL_LINE_STRIP)

	gl.glVertex3f(chart:chart(cartesianToLatLonWGS84(pos:unpack())))
	for i=1,numSteps do

		q:rotate(pos, pos)
--print(pos, chart:chart(cartesianToLatLonWGS84(pos:unpack())))
		gl.glVertex3f(chart:chart(cartesianToLatLonWGS84(pos:unpack())))
	end

	--[[
	ok how to do this?  geodesics and connections?  rotate spherical coordinates for now?
	so
	1) get the basis at this point, with one basis vector pointing to true (rotation axis) north, and another at 90' of it
	2) project our magnetic north heading onto this axis
	3) that will give us a geodesic direction
	4) cross that with our cartesian position to get our axis
	5) rotate to get our great arc
	--]]

	gl.glEnd()
end

local function drawVectorField(chart, app)

	local arrow = {
		-- [[ arrow for real
		{-.5, 0.},
		{.5, 0.},
		{.2, .3},
		{.5, 0.},
		{.2, -.3},
		{.5, 0.},
		--]]
		--[[ debug basis display
		{0,0,0},
		{1,0,0},
		{0,0,0},
		{0,1,0},
		{0,0,0},
		{0,0,1},
		--]]
	}

	local height = 0
	local londim = 60
	local latdim = 30
	-- [[
	local scale = guivars.arrowScale  / (BStat.mag.max * latdim)
	--]]
	--[[ debug basis display
	local scale = .5 / 30
	--]]

	local shader = app.vectorFieldShader
	shader:use()
	shader:setUniforms{
		dt = guivars.fieldDT,
		mvProjMat = app.view.mvProjMat.ptr,
		arrowScale = scale,
		fieldZScale = guivars.fieldZScale,
		weight_Equirectangular = guivars.geomIndex == chartIndexForName.Equirectangular and 1 or 0,
		weight_Azimuthal_equidistant = guivars.geomIndex == chartIndexForName['Azimuthal equidistant'] and 1 or 0,
		weight_Mollweide = guivars.geomIndex == chartIndexForName.Mollweide and 1 or 0,
		weight_WGS84 = guivars.geomIndex == chartIndexForName.WGS84 and 1 or 0,
		chartIs3D = chart.is3D or false,
	}
	chart.basisTex:bind()
	--[[ B coeffs using textures
	BTex:bind()
	--]]

	gl.glBegin(gl.GL_LINES)

	for j=0,latdim-1 do
		local v = (j + .5) / latdim
		local lat = (v * 2 - 1) * 90
		local phi = math.rad(lat)

		--[[ equal res at all latitudes
		local thislondim = londim
		--]]
		-- [[ proportional to the latitude arclen
		local thislondim = math.max(1, math.ceil(londim * math.cos(phi)))
		--]]
		for i=0,thislondim-1 do
			local u = (i + .5) / thislondim
			local lon = (u * 2 - 1) * 180
			local lambda = math.rad(lon)

			gl.glVertexAttrib2f(shader.attrs.texcoord.loc, u, v)

			-- [[
			local ex, ey, ez = chart:basis(phi, lambda, height)
			gl.glVertexAttrib3f(shader.attrs.ex.loc, ex:unpack())
			gl.glVertexAttrib3f(shader.attrs.ey.loc, ey:unpack())
			gl.glVertexAttrib3f(shader.attrs.ez.loc, ez:unpack())
			--]]

			-- TODO instanced geometry?  any faster?
			for _,q in ipairs(arrow) do
				gl.glVertex2f(q[1], q[2])
			end
		end
	end
	gl.glEnd()

	GLTex2D:unbind()
	shader:useNone()
end


function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local shader = assert(overlays[tonumber(guivars.overlayIndex)+1]).shader
	local gradtex = assert(gradients[tonumber(guivars.gradientIndex)+1]).tex
	local chart = assert(charts[tonumber(guivars.geomIndex)])

	shader:use()

	if shader.uniforms.dt then
		gl.glUniform1f(shader.uniforms.dt.loc, guivars.fieldDT)
	end

	earthtex:bind(0)
	gradtex:bind(1)
	BTex:bind(2)
	B2Tex:bind(3)

	gl.glUniform1f(shader.uniforms.alpha.loc, 1)
	gl.glUniform1f(shader.uniforms.gradScale.loc, guivars.gradScale)

	shader:useNone()

-- TODO more samples
--[[
	drawReading{
		chart = chart,
		lat = {42, 52, 45.021},
		lon = {74, 34, 16.004},
		heading = 5,
	}

	drawReading{
		chart = chart,
		lat = {33, 59, 38},		-- lat
		lon = {-80, -27, -56},	-- lon
		heading = -.5,
	}
--]]

	if guivars.doDrawVectorField then
		drawVectorField(chart, self)
	end

	shader:use()

	gl.glUniform1f(shader.uniforms.alpha.loc, guivars.drawAlpha)

	gl.glCullFace(gl.GL_FRONT)
	chart:draw(self, shader, gradtex)
	gl.glCullFace(gl.GL_BACK)
	chart:draw(self, shader, gradtex)

	B2Tex:unbind(3)
	BTex:unbind(2)
	gradtex:unbind(1)
	earthtex:unbind(0)

	shader:useNone()

	glreport'here'

	--render gui
	App.super.update(self, ...)
end

function App:updateGUI()

	local thisTime = os.time()
	if thisTime ~= self.lastTime then
		if self.lastTime then
			self.fps = self.frames / (thisTime - self.lastTime)
			self.frames = 0
		end
		self.lastTime = thisTime
	end
	self.frames = (self.frames or 0) + 1

	ig.igText('fps: '..tostring(self.fps))


	ig.luatableCheckbox('ortho', self.view, 'ortho')
	ig.luatableCheckbox('draw vector field', guivars, 'doDrawVectorField')

	ig.igText'chart'
	for i,chart in ipairs(charts) do
		ig.luatableRadioButton(chart.name, guivars, 'geomIndex', i)
	end

	ig.igSeparator()

	local chart = assert(charts[tonumber(guivars.geomIndex)])
	if chart.updateGUI then
		chart:updateGUI()
	end

	ig.igSeparator()

	ig.igText'overlay'
	for i,overlay in ipairs(overlays) do
		ig.luatableRadioButton(overlay.name, guivars, 'overlayIndex', i-1)
	end

	ig.igSeparator()

	ig.igText'gradient'
	for i,grad in ipairs(gradients) do
		ig.luatableRadioButton(grad.name, guivars, 'gradientIndex', i-1)
	end
	ig.luatableInputFloat('gradScale', guivars, 'gradScale')

	ig.luatableInputFloat('alpha', guivars, 'drawAlpha')
	ig.luatableInputFloat('field z', guivars, 'fieldZScale')
	ig.luatableInputFloat('field size', guivars, 'arrowScale')

	-- how linear are the g and h coeffs?
	-- can I just factor out the dt?
	--ig.luatableInputFloat('time from '..wmm.epoch, guivars, 'fieldDT')
	ig.luatableSliderFloat('years from '..tostring(wmm.epoch), guivars, 'fieldDT', -50, 50)
end

return App():run()

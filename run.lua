#!/usr/bin/env luajit
local cmdline = require 'ext.cmdline'(...)
local assertindex = require 'ext.assert'.index
local math = require 'ext.math'	-- isfinite
local timer = require 'ext.timer'
local gl = require 'gl.setup'(cmdline.gl or 'OpenGL')
local ffi = require 'ffi'
local vec2f = require 'vec-ffi.vec2f'
local vec3f = require 'vec-ffi.vec3f'
local vec4f = require 'vec-ffi.vec4f'
local quatf = require 'vec-ffi.quatf'
local vector = require 'ffi.cpp.vector-lua'
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

-- ported from WMM2020 GeomagnetismLibrary.c
-- phi = radians
-- lambda = radians
-- height = meters
local function calcB(phi, lambda, height)
	height = height * 1e-3	-- m to km

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

-- ported from WMM2020 GeomagnetismLibrary.c
-- expects xyz in cartesian units earth-semimajor-axis
-- output is in (radians, radians, km)
-- TODO just use charts.WGS84:chartInv(x,y,z) ?  or use this there?
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
	return phi, lambda, height * 1e+3		-- km back to m
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
	overlayIndex = 2,
	gradientIndex = 1,

	drawAlpha = 1,
	doDrawArrowField = true,
	doDrawFieldLines = true,

	fieldDT = 0,

	arrowHeight = 0,
	arrowLatDim = 30,
	arrowLonDim = 60,
	arrowEvenDistribute = true,
	arrowScale = 5,
	-- set to 0 to flatten vector field against surface
	arrowZScale = 1,
	arrowFieldNormalize = false,

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
		return oldBasisFunc(self, math.deg(phi), math.deg(lambda), height)
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
	function c:draw(app, shader, gradTex)
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
				texs = {earthtex, gradTex, BTex, B2Tex},
			}
		end

		self.sceneobj.program = shader
		self.sceneobj.texs[2] = gradTex
		self.sceneobj.uniforms.mvProjMat = app.view.mvProjMat.ptr
		self.sceneobj.uniforms.dt = guivars.fieldDT
		self.sceneobj.uniforms.alpha = guivars.drawAlpha
		self.sceneobj.uniforms.gradScale = guivars.gradScale
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

	local londim = 1440
	local latdim = 720

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

	local BData = ffi.new('vec4f_t[?]', londim * latdim)
	local B2Data = ffi.new('vec4f_t[?]', londim * latdim)
	timer('generating B field', function()
		-- can be arbitrary
		-- but the WMM model is for 15 mins, so [360,180] x4

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
			program = {
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
			},
			geometry = self.quadGeom,
		}

		-- BData / B2Data is only used for stat computation
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
	end)

--[[
		for j=0,latdim-1 do
			local v = (j + .5) / latdim
			local lat = (v * 2 - 1) * 90
			local phi = math.rad(lat)
			for i=0,londim-1 do
				local u = (i + .5) / londim
				local lon = (u * 2 - 1) * 180
				local lambda = math.rad(lon)
				print(BData[i+londim*j], calcB(phi, lambda, 0))
			end
		end
--]]

-- [=[ hmm, better way than copy paste?

	timer('generating div B and curl B', function()
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
			program = {
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
vec3 dheight = vec3(0., 0., 1000.);

in vec2 texcoordv;
out vec4 fragColor;

void main() {
	float phi = (texcoordv.y - .5) * M_PI;			//[-pi/2, pi/2]
	float lambda = (texcoordv.x - .5) * 2. * M_PI;	//[-pi, pi]

	vec3 plh = vec3(phi, lambda, 0.);

	// TODO units anyone?
	vec3 dphi_B = (calcB(plh + dphi) - calcB(plh - dphi)) / dphi.x / (wgs84_a * 1e+3 * cos(plh.x));
	vec3 dlambda_B = (calcB(plh + dlambda) - calcB(plh - dlambda)) / dlambda.y / (wgs84_a * 1e+3);
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
]],				{
					latdim = latdim,
					londim = londim,
					clnumber = clnumber,
				}),
				uniforms = {
					dt = 0,
				},
			},
			geometry = self.quadGeom,
		}

		-- only used for stat calc
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
	end)

	timer('generating stats', function()
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
	end)

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

	self.view.ortho = true
	self.view.pos.z = 2
	self.view.orthoSize = 2
	gl.glClearColor(0,0,0,0)

--[=[
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
--]=]

	timer('building arrow field', function()
		local arrowFieldShader = GLProgram{
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

layout(location=0) in vec2 vertex;

uniform mat4 mvProjMat;

uniform float BMagMax;

uniform float arrowHeight;
uniform float arrowLatDim;
uniform float arrowLonDim;
uniform float arrowEvenDistribute;
uniform float arrowScale;
uniform float arrowZScale;
uniform bool arrowFieldNormalize;

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

uniform int arrowTexSize;
uniform sampler2D arrowCoordTex;

void main() {
	vec2 arrowTexTC;
	int iarrowTexTCX = gl_InstanceID % arrowTexSize;
	int iarrowTexTCY = (gl_InstanceID - iarrowTexTCX) / arrowTexSize;
	arrowTexTC.x = float(iarrowTexTCX);	// integer in x direction
	arrowTexTC.y = float(iarrowTexTCY);	// integer in y direction
	arrowTexTC += .5;
	arrowTexTC /= float(arrowTexSize);
	vec2 texcoord = texture(arrowCoordTex, arrowTexTC).xy;

	vec3 latLonHeight = vec3(
		(texcoord.y - .5) * 180.,	// lat in deg
		(texcoord.x - .5) * 360., 	// lon in deg
		arrowHeight);				// height in m

	// expect vertex xyz to be lat lon height
	// lat and lon is in degrees
	// height is in meters
	vec3 pos = 0.
<? for _,name in ipairs(chartCNames) do
?>		+ weight_<?=name?> * chart_<?=name?>(latLonHeight)
<? end
?>	;
	//from meters to normalized coordinates
	pos /= WGS84_a;

	if (chartIs3D) {
		pos = vec3(pos.z, pos.x, pos.y);
	}

	vec3 BCoeffs = calcB(vec3(
		(texcoord.y - .5) * M_PI,			//phi
		(texcoord.x - .5) * 2. * M_PI,		//lambda
		0.
	));

	mat3 e = mat3(vec3(0.), vec3(0.), vec3(0.))
<? for _,name in ipairs(chartCNames) do
?>		+ weight_<?=name?> * chart_<?=name?>_basis(latLonHeight)
<? end
?>	;

#if 0	// see the comments in Chart:getGLSLFunc3D for why I'm permuting basis coords there and not here ...
	if (chartIs3D) {
		e[0] = xformZBackToZUp(e[0]);
		e[1] = xformZBackToZUp(e[1]);
		e[2] = xformZBackToZUp(e[2]);
	}
#endif

	// do this before normalizing
	BCoeffs.z *= arrowZScale;

	//should this normalize after applying scales/basis or before?
	// before, so the scales can do something
	if (arrowFieldNormalize) {
		BCoeffs = normalize(BCoeffs);
	} else {
		BCoeffs = BCoeffs / BMagMax;
	}

	vec3 B = (arrowScale / arrowLatDim) * (e * BCoeffs);	// e_i B^i

	//vec3 Bey = perpTo(B);			// cross with furthest axis
	//vec3 Bey = vec3(-B.y, B.x, 0.);	// cross with z axis
	vec3 Bey = (arrowScale / arrowLatDim) * (e * vec3(-BCoeffs.y, BCoeffs.x, BCoeffs.z));

	gl_Position = mvProjMat * vec4(
		pos + vertex.x * B
			+ vertex.y * Bey,
		1.);
}
]],
			{
				clnumber = clnumber,
				chartCNames = chartCNames,
			}),
			fragmentCode = [[
out vec4 fragColor;
void main() {
	fragColor = vec4(1., 1., 1., 1.);
}
]],
		}:useNone()

		self.arrowFieldSceneObj = GLSceneObject{
			program = arrowFieldShader,
			vertexes = {
				data = {
					-.5, 0,
					.5, 0,
					.2, .3,
					.5, 0,
					.2, -.3,
					.5, 0,
				},
				dim = 2,
			},
			geometry = {
				mode = gl.GL_LINES,
			},
			uniforms = {
				arrowTexSize = arrowTexSize,
			},
		}

		self:updateArrowTex()	-- do this after self.arrowFieldSceneObj is made
	end)

	timer('building field lines', function()
		local fieldLineShader = GLProgram{
			version = 'latest',
			precision = 'best',
			vertexCode = chartCode..template([[
<? for _,name in ipairs(chartCNames) do
?>uniform float weight_<?=name?>;
<? end
?>
uniform bool chartIs3D;

uniform vec2 BMag;	// (min, max)

//vertex = (phi (rad), lambda (rad), height (meters), |B|)
#if 1 // vertexes from attrib
layout(location=0) in vec4 vertex;
#endif

#if 0 // vertexes from texture
uniform int fieldLineVtxsTexWidth;	// equals iterations in a line strip
uniform int fieldLineVtxsTexHeight;	// equals # of start positions
uniform sampler2D fieldLineVtxsTex;
#endif

out float BFrac;

uniform mat4 mvProjMat;
void main() {
#if 0 // vertexes from texture
	vec2 fieldLineVtxsTC;
	int indexInLine = gl_VertexID % fieldLineVtxsTexWidth;
	fieldLineVtxsTC.x = (float(indexInLine) + .5) / float(fieldLineVtxsTexWidth);
	fieldLineVtxsTC.y = (float((gl_VertexID - indexInLine) / fieldLineVtxsTexWidth) + .5) / float(fieldLineVtxsTexHeight);
	vec4 vertex = texture(fieldLineVtxsTex, fieldLineVtxsTC);
#endif

	vec3 latLonHeight = vec3(
		vertex.xy * (180. / M_PI),
		vertex.z);

	BFrac = (vertex.w - BMag.x) / (BMag.y - BMag.x);

	vec3 pos = 0.
<? for _,name in ipairs(chartCNames) do
?>		+ weight_<?=name?> * chart_<?=name?>(latLonHeight)
<? end
?>	;
	//from meters to normalized coordinates
	pos /= WGS84_a;

	if (chartIs3D) {
		pos = vec3(pos.z, pos.x, pos.y);
	}

	gl_Position = mvProjMat * vec4(pos, 1.);
}
]],			{
				chartCNames = chartCNames,
			}),
			fragmentCode = [[
in float BFrac;
out vec4 fragColor;
uniform sampler2D gradTex;
void main() {
	fragColor = texture(gradTex, vec2(BFrac, .5));
}
]],
			{
				gradTex = 0,
				fieldLineVtxsTex = 1,
			},
		}:useNone()

		local londim = 30	-- field line lon dim
		local latdim = 15	-- field line lat dim
		local evenDistribute = false	-- field line even distribute
		local height0 = 0	-- field line initial height

--[[
how to build the field lines on the GPU ...
1) pick our seed locations / fill them in on the top row of a float tex
2) run a 1 x n shader / fbo on each successive rows, integrate as we go ...
3) copy the fbo tex to an array buffer, and use it as a bunch of line strip vertexes
--]]

		local startCoords = vector'vec4f_t'
		-- try to draw magnetic field lines ...
		for j=0,latdim-1 do
			local v = (j + .5) / latdim
			local lat = (v * 2 - 1) * 90
			local phi = math.rad(lat)

			local thislondim = evenDistribute
				and math.max(1, math.ceil(londim * math.cos(phi)))
				or londim
			for i=0,thislondim-1 do
				local u = (i + .5) / thislondim
				local lon = (u * 2 - 1) * 180
				local lambda = math.rad(lon)

				-- TODO Bmag on the GPU would be nice ...
				local Bx, By, Bz = calcB(phi, lambda, height0)
				local Bmag = math.sqrt(Bx^2 + By^2 + Bz^2)

				for sign=1,2 do
					-- assert here that every other startCoords is going to integrate + - + - ...
					startCoords:emplace_back()[0]:set(phi, lambda, height0, Bmag)
				end
			end
		end
print('# field line start coords', #startCoords)
		
		local dparam = 1e-2
		local maxiter = 400

-- [==[ new way, integrate on GPU
		-- store current state here
		-- I couldve used a pingpong but why waste the memory
		local fieldLinePosTex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,	-- phi, lambda, height, |B|
			width = 1,			-- 1D along height for blitting's sake 
			height = #startCoords,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_NEAREST,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
			data = startCoords.v,
		}:unbind()

		-- use a fbo to integrate / generate our field line positions
		-- seed on the y axis so that successive lines form along the x-axis i.e. along contiguous memory to be later passed off to the draw calls
		-- (use a PBO?)
		self.fieldLineVtxsTex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,	-- phi, lambda, height, |B|
			width = maxiter,
			height = #startCoords,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_NEAREST,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}:subimage{	-- upload vertical strip
			width = 1,
			height = #startCoords,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			data = startCoords.v,
		}:unbind()

		local calcFieldLineVtxsTexSceneObj = GLSceneObject{
			program = {
				version = 'latest',
				precision = 'best',
				vertexCode = quadGeomVertexCode,
				fragmentCode = chartCode
..[[
uniform float dt;
]]..calc_b_shader..[[

// ported from WMM2020 GeomagnetismLibrary.c
// expects xyz in cartesian units earth-semimajor-axis
// output is in (radians, radians, km)
// TODO just use charts.WGS84:chartInv(x,y,z) ?  or use this there?
vec3 cartesianToLatLonWGS84(vec3 pos) {
	// lowercase wgs84_a is in km, uppercase WGS84_a is in m ... i know i need to fix this ...
	pos *= wgs84_a;	// convert from semimajor-axis units to km
	float modified_b = pos.z < 0 ? -wgs84_b : wgs84_b;
	float r = length(pos.xy);

	float e = (modified_b * pos.z - (wgs84_a * wgs84_a - modified_b * modified_b)) / (wgs84_a * r);
	float f = (modified_b * pos.z + (wgs84_a * wgs84_a - modified_b * modified_b)) / (wgs84_a * r);
	float p = (4. / 3.) * (e * f + 1.);
	float q = 2. * (e * e - f * f);
	float d = p * p * p + q * q;

	float v;
	if  (d >= 0) {
		v = pow(sqrt(d) - q, 1./3.) - pow(sqrt(d) + q, 1./3.);
	} else {
		v = 2. * sqrt(-p) * cos(acos(q / (p * sqrt(-p))) / 3.);
	}

	if (v * v < abs(p)) {
		v = -(v * v * v + 2. * q) / (3. * p);
	}

	float g = (sqrt(e * e + v) + e) / 2.;
	float t = sqrt(g * g + (f - v * g) / (2. * g - e)) - g;

	float rlat = atan((wgs84_a * (1. - t * t)) / (2. * modified_b * t));
	float phi = rlat;

	float height = (r - wgs84_a * t) * cos(rlat) + (pos.z - modified_b) * sin(rlat);
	float zlong = atan(pos.y, pos.x);
	if (zlong < 0.) {
		zlong += 2.*M_PI;
	}
	float lambda = zlong;
	lambda += M_PI;
	lambda = mod(lambda, 2. * M_PI);
	lambda -= M_PI;
	return vec3(phi, lambda, height * 1e+3);		// km back to m
}

in vec2 texcoordv;
out vec4 fragColor;
uniform float dparam;
uniform float texHeight;
uniform sampler2D fieldLinePosTex;

void main() {
	// sign is -1 for even rows, 1 for odd rows
	float sign = floor(mod(texcoordv.y * texHeight, 2.)) * 2. - 1.;

	// texcoordv.x is the iteration, texcoordv.y is the coord index
	// meanwhile fieldLinePosTex.x is the coord index, and it has no .y
	vec4 phiLambdaHeightBMag = texture(fieldLinePosTex, vec2(.5, texcoordv.y));
	vec3 latLonHeight = vec3(
		phiLambdaHeightBMag.xy * 180. / M_PI, 
		phiLambdaHeightBMag.z);

	// calculate cartesian position, returns in meters, divide by WGS84_a (in meters) to get cartesian coords in units of earth semimajor axis
	vec3 cartesianPos = chart_WGS84(latLonHeight) / WGS84_a;
cartesianPos = vec3(cartesianPos.z, cartesianPos.x, cartesianPos.y);	// undo xformZBackToZUp

// but don't undo the basis xform?  why ....
	mat3 e = chart_WGS84_basis(latLonHeight);

	vec3 B = calcB(phiLambdaHeightBMag.xyz);	// input: phi, lambda, height
	
	// should match phiLambdaHeightBMag.w which is the last bmag
	// maybe I could save calcs by storing a second tex of the bvec alongside the plh vec ...
	//vec3 BMag = length(B);

	vec3 dv = e * B;//B * e;	// e * B or B * e ?
	vec3 n = normalize(dv) * sign * dparam;
	cartesianPos += n;
	//if (length(cartesianPos) < .02) then stop integrating or something meh

	phiLambdaHeightBMag.xyz = cartesianToLatLonWGS84(cartesianPos);
	//and recompute the new coord's B and store the new Bmag ... wasted calcs ... I could avoid by storing both phi,lambda,hegiht and storing Bvec as we integrate,  but that's one extra texture ...
	B = calcB(phiLambdaHeightBMag.xyz);
	phiLambdaHeightBMag.w = length(B);
	fragColor = phiLambdaHeightBMag;
}
]],
				uniforms = {
					dt = 0,
					dparam = dparam,
					texHeight = #startCoords,
					fieldLinePosTex = 0,
				},
			},
			geometry = self.quadGeom,
			texs = {fieldLinePosTex},
		}

		local stateFBO = GLFBO()
			:setColorAttachmentTex2D(fieldLinePosTex.id)
		local res, err = fbo.check()
		if not res then print(err) end
		stateFBO:unbind()

		timer('integrating on GPU', function()
			-- now iteratively draw to FBO next strip over, reading from the current state strip as we go
			-- (TODO store a current-state tex separately?)
			for i=1,maxiter-1 do
				-- can I keep fbo bound as GL_FRAMEBUFFER and as GL_DRAW_FRAMEBUFFER?
				fbo:bind()
					:setColorAttachmentTex2D(self.fieldLineVtxsTex.id)
				local res, err = fbo.check()
				if not res then print(err) end
				-- draw from fieldLinePosTex into fieldLineVtxsTex
				gl.glViewport(i, 0, 1, self.fieldLineVtxsTex.height)
				calcFieldLineVtxsTexSceneObj.uniforms.mvProjMat = self.unitProjMat.ptr
				calcFieldLineVtxsTexSceneObj:draw()
				fbo:unbind()

				-- read back the i'th column into the fieldLinePosTex texture
				-- draw from fieldLineVtxsTex into fieldLinePosTex
				gl.glBindFramebuffer(gl.GL_READ_FRAMEBUFFER, fbo.id)
				gl.glBindFramebuffer(gl.GL_DRAW_FRAMEBUFFER, stateFBO.id)
				-- readpixels wrt viewport or framebuffer origin?
				-- means use pixel buffers?  or use pixel buffers once i get around to writing this to vertex buffers for the geom?
				--gl.glReadPixels(i, 0, 1, self.fieldLineVtxsTex.height, gl.GL_RGBA, gl.GL_FLOAT, ffi.cast('void*', 0));
				gl.glBlitFramebuffer(
					i,		-- srcX0
					0,		-- srcY0
					i+1,	-- srcX1
					self.fieldLineVtxsTex.height,	-- srcY1
					0,		-- dstX0
					0,		-- dstY0
					1,		-- dstX1
					self.fieldLineVtxsTex.height,	-- dstY1
					gl.GL_COLOR_BUFFER_BIT,
					gl.GL_NEAREST)
				gl.glBindFramebuffer(gl.GL_READ_FRAMEBUFFER, 0)
				gl.glBindFramebuffer(gl.GL_DRAW_FRAMEBUFFER, 0)
			end
		end)
--]==]

--[==[ upload with cpu copy
		local vertexes = table()
		local geometries = table()

		local fieldLineVtxsData = ffi.new('vec4f_t[?]', maxiter * #startCoords)
		self.fieldLineVtxsTex:toCPU(fieldLineVtxsData)

		for startIndex=0,#startCoords-1 do
			local indexStart = #vertexes / 4

			for k=0,maxiter-1 do
				local src = fieldLineVtxsData[k + maxiter * startIndex]
				vertexes:insert(src.x)
				vertexes:insert(src.y)
				vertexes:insert(src.z)
				vertexes:insert(src.w)
			end
			local indexEnd = #vertexes / 4

			-- giving each geometry a different .vertexes array would remove the need for indicies
			-- and on old webgl that would mean getting rid of a 65536 max vertex limit
			-- but that'd mean changing sceneobject to rebind the 'vertexes' attribute each time the geometry changes
			-- which maybe I should do since the .geometries field is new
			-- and this design could still benefit from glMultiDrawArrays ...
			-- ... or not, that'd have to be an explicit separate case, since it requires a single vertexes to be bound.
			geometries:insert{
				mode = gl.GL_LINE_STRIP,
				count = indexEnd - indexStart,
				offset = indexStart,
			}
		end

		self.fieldLineSceneObj = GLSceneObject{
			program = fieldLineShader,
			vertexes = {
				data = vertexes,
				dim = 4,
			},
			geometries = geometries,
		}
--]==]
-- [==[ set vertexes via gpu
	
		-- [===[ can't seem to get pbos to work ....
		-- library design fallacy, unlike textures, the same gl buffer can be bound to various targets
		self.fieldLineVertexBuf = require 'gl.pixelpackbuffer'{
			-- [[
			size = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height * 4 * 4,	-- sizeof(float) * 4 components
			type = gl.GL_FLOAT,
			count = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height,
			dim = 4,
			mode = gl.GL_DYNAMIC_DRAW,--gl.GL_STATIC_DRAW,
			--]]
		}
		--[[
		gl.glBufferData(gl.GL_PIXEL_PACK_BUFFER, self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height * 4 * 4, nil, gl.GL_DYNAMIC_DRAW)
		--self.fieldLineVertexBuf:unbind()
		-- assign afterwards or the class will allocate filler
		self.fieldLineVertexBuf.size = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height * 4 * 4
		self.fieldLineVertexBuf.type = gl.GL_FLOAT
		self.fieldLineVertexBuf.count = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height
		self.fieldLineVertexBuf.dim = 4
		self.fieldLineVertexBuf.mode = gl.GL_DYNAMIC_DRAW
		--self.fieldLineVertexBuf.mode = gl.GL_STATIC_DRAW
		--]]

glreport'here'
		fbo:bind()
glreport'here'
		fbo:setColorAttachmentTex2D(self.fieldLineVtxsTex.id)
glreport'here'
		local res, err = fbo.check()
		if not res then print(err) end
glreport'here'

glreport'here'
		gl.glViewport(0, 0, self.fieldLineVtxsTex.width, self.fieldLineVtxsTex.height)
		-- read from bound FBO (attached to tex) to PBO (attached to buffer)
		gl.glReadPixels(0, 0, self.fieldLineVtxsTex.width, self.fieldLineVtxsTex.height, gl.GL_RGBA, gl.GL_FLOAT, ffi.cast('void*', 0))
glreport'here'
		fbo:unbind()
glreport'here'
		self.fieldLineVertexBuf:unbind()
glreport'here'

		-- and now the PBO buffer can be treated as an array buffer...so...
		--self.fieldLineVertexBuf.target = gl.GL_ARRAY_BUFFER
		-- don't just change the target, change the class, so the sceneobject vertexes class-detect will auto-assign it to .buffer
		setmetatable(self.fieldLineVertexBuf, require 'gl.arraybuffer')
		--self.fieldLineVertexBuf:bind():unbind()
		--]===]
		--[===[ how about unpack buffer + texsubimage
		self.fieldLineVertexBuf = require 'gl.pixelunpackbuffer'{
			size = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height * 4 * 4,	-- sizeof(float) * 4 components
			type = gl.GL_FLOAT,
			count = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height,
			dim = 4,
			mode = gl.GL_DYNAMIC_DRAW,--gl.GL_STATIC_DRAW,
		}
		self.fieldLineVtxsTex
			:bind()
			:subimage{x=0, y=0, width=self.fieldLineVtxsTex.width, height=self.fieldLineVtxsTex.height, format=gl.GL_RGBA, type=gl.GL_FLOAT, data=ffi.cast('void*', 0)}
			:unbind()
		self.fieldLineVertexBuf:unbind()
		--self.fieldLineVertexBuf.target = gl.GL_ARRAY_BUFFER
		-- don't just change the target, change the class, so the sceneobject vertexes class-detect will auto-assign it to .buffer
		setmetatable(self.fieldLineVertexBuf, require 'gl.arraybuffer')
		--]===]

		local geometries = table()
		for i=0,#startCoords-1 do
			geometries:insert{
				mode = gl.GL_LINE_STRIP,
				count = maxiter,
				offset = i * maxiter,
			}
		end
		
		self.fieldLineSceneObj = GLSceneObject{
			program = fieldLineShader,
			vertexes = self.fieldLineVertexBuf,
			geometries = geometries,
			--[[
			uniforms = {
				fieldLineVtxsTexWidth = self.fieldLineVtxsTex.width,
				fieldLineVtxsTexHeight = self.fieldLineVtxsTex.height,
			},
			--]]
		}
glreport'here'
--]==]
	end)
	
	--gl.glEnable(gl.GL_CULL_FACE)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
end

function App:updateArrowTex()
	-- whenever lat/lon are resized, just rebuild this
	-- technically you could even do it on the GPU as a GPU pass ...
	-- just 1) unravel the texcoord to 1D index
	-- then 2) ... well the inverse map of the sum of ceil(londim * cos(phi)) is a bit tricky ...
	local arrowCoords = vector'vec2f_t'
	local londim = guivars.arrowLonDim
	local latdim = guivars.arrowLatDim
	for j=0,latdim-1 do
		local v = (j + .5) / latdim
		local lat = (v * 2 - 1) * 90
		local phi = math.rad(lat)

		local thislondim = guivars.arrowEvenDistribute
			and math.max(1, math.ceil(londim * math.cos(phi)))
			or londim
		for i=0,thislondim-1 do
			local u = (i + .5) / thislondim
			local lon = (u * 2 - 1) * 180
			local lambda = math.rad(lon)
			arrowCoords:emplace_back()[0]:set(u, v)
		end
	end
	local arrowTexSize = math.ceil(math.sqrt(#arrowCoords))
	arrowCoords:reserve(arrowTexSize * arrowTexSize)

	local arrowCoordTex = GLTex2D{
		internalFormat = gl.GL_RG32F,
		width = arrowTexSize,
		height = arrowTexSize,
		format = gl.GL_RG,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_NEAREST,
		wrap = {
			s = gl.GL_REPEAT,
			t = gl.GL_REPEAT,
		},
		data = arrowCoords.v,
	}:unbind()

	self.arrowFieldSceneObj.geometry.instanceCount = #arrowCoords
	self.arrowFieldSceneObj.uniforms.arrowTexSize = arrowTexSize
	self.arrowFieldSceneObj.texs[1] = arrowCoordTex
end

local function degmintofrac(deg, min, sec)
	return deg + (1/60) * (min + (1/60) * sec)
end

--[=[
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
	local pos = vec3f(charts.WGS84:chart(math.deg(phi), math.deg(lambda), height)) / allCharts.WGS84_a
	local ex, ey, ez = charts.WGS84:basis(math.deg(phi), math.deg(lambda), .1)
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
--]=]

function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local shader = assert(overlays[tonumber(guivars.overlayIndex)]).shader
	local gradTex = assert(gradients[tonumber(guivars.gradientIndex)]).tex
	local chart = assert(charts[tonumber(guivars.geomIndex)])

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

	if guivars.doDrawArrowField then
		self.arrowFieldSceneObj.uniforms.dt = guivars.fieldDT
		self.arrowFieldSceneObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
		self.arrowFieldSceneObj.uniforms.arrowZScale = guivars.arrowZScale
		self.arrowFieldSceneObj.uniforms.arrowFieldNormalize = guivars.arrowFieldNormalize
		self.arrowFieldSceneObj.uniforms.arrowHeight = guivars.arrowHeight
		self.arrowFieldSceneObj.uniforms.arrowLatDim = guivars.arrowLatDim
		self.arrowFieldSceneObj.uniforms.arrowLonDim = guivars.arrowLonDim
		self.arrowFieldSceneObj.uniforms.arrowEvenDistribute = guivars.arrowEvenDistribute
		self.arrowFieldSceneObj.uniforms.BMagMax = BStat.mag.max
		self.arrowFieldSceneObj.uniforms.arrowScale = guivars.arrowScale

		self.arrowFieldSceneObj.uniforms.weight_Equirectangular = guivars.geomIndex == chartIndexForName.Equirectangular and 1 or 0
		self.arrowFieldSceneObj.uniforms.weight_Azimuthal_equidistant = guivars.geomIndex == chartIndexForName['Azimuthal equidistant'] and 1 or 0
		self.arrowFieldSceneObj.uniforms.weight_Mollweide = guivars.geomIndex == chartIndexForName.Mollweide and 1 or 0
		self.arrowFieldSceneObj.uniforms.weight_WGS84 = guivars.geomIndex == chartIndexForName.WGS84 and 1 or 0
		self.arrowFieldSceneObj.uniforms.chartIs3D = chart.is3D or false

		self.arrowFieldSceneObj:draw()
	end
	if guivars.doDrawFieldLines then
		self.fieldLineSceneObj.texs[1] = gradTex
		self.fieldLineSceneObj.texs[2] = self.fieldLineVtxsTex
		self.fieldLineSceneObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
		self.fieldLineSceneObj.uniforms.BMag = {BStat.mag.min, BStat.mag.max}

		-- TODO if we're using weights then why do we need the 'chartIs3D' flag?
		self.fieldLineSceneObj.uniforms.weight_Equirectangular = guivars.geomIndex == chartIndexForName.Equirectangular and 1 or 0
		self.fieldLineSceneObj.uniforms.weight_Azimuthal_equidistant = guivars.geomIndex == chartIndexForName['Azimuthal equidistant'] and 1 or 0
		self.fieldLineSceneObj.uniforms.weight_Mollweide = guivars.geomIndex == chartIndexForName.Mollweide and 1 or 0
		self.fieldLineSceneObj.uniforms.weight_WGS84 = guivars.geomIndex == chartIndexForName.WGS84 and 1 or 0
		self.fieldLineSceneObj.uniforms.chartIs3D = chart.is3D or false

		self.fieldLineSceneObj:draw()
	end

	-- why cull each side separately?
	--gl.glCullFace(gl.GL_FRONT)
	chart:draw(self, shader, gradTex)
	--gl.glCullFace(gl.GL_BACK)
	--chart:draw(self, shader, gradTex)

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
	ig.luatableCheckbox('draw arrow field', guivars, 'doDrawArrowField')
	ig.luatableCheckbox('draw field lines', guivars, 'doDrawFieldLines')

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
		ig.luatableRadioButton(overlay.name, guivars, 'overlayIndex', i)
	end

	ig.igSeparator()

	ig.igText'gradient'
	for i,grad in ipairs(gradients) do
		ig.luatableRadioButton(grad.name, guivars, 'gradientIndex', i)
	end
	ig.luatableInputFloat('gradScale', guivars, 'gradScale')

	ig.luatableInputFloat('alpha', guivars, 'drawAlpha')
	ig.luatableInputFloat('arrow altitude', guivars, 'arrowHeight')
	ig.luatableInputFloat('arrow scale', guivars, 'arrowScale')
	ig.luatableInputFloat('arrow z scale', guivars, 'arrowZScale')
	ig.luatableCheckbox('arrow field normalize', guivars, 'arrowFieldNormalize')

	if ig.luatableInputInt('arrow lat dim', guivars, 'arrowLatDim')
	or ig.luatableInputInt('arrow lon dim', guivars, 'arrowLonDim')
	or ig.luatableCheckbox('arrows even distribute', guivars, 'arrowEvenDistribute')
	then
		self:updateArrowTex()
	end

	-- how linear are the g and h coeffs?
	-- can I just factor out the dt?
	--ig.luatableInputFloat('time from '..wmm.epoch, guivars, 'fieldDT')
	ig.luatableSliderFloat('years from '..tostring(wmm.epoch), guivars, 'fieldDT', -50, 50)
end

return App():run()

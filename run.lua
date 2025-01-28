#!/usr/bin/env luajit
local cmdline = require 'ext.cmdline'(...)
local assert = require 'ext.assert'
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
local range = require 'ext.range'
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

local wmmfn = cmdline.cof or 'wmm2020/wmm.cof'

-- load wmm
local lines = string.split(string.trim((assert(path(wmmfn):read()))), '\n')
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

local earthtex

local App = require 'imguiapp.withorbit'()
App.title = 'EM field'


local guivars = {
	geomIndex = 1,
	overlayIndex = 2,
	gradientIndex = 1,

	drawAlpha = 1,
	doDrawArrowField = true,
	doDrawFieldLines = true,

	gradScale = 1,

	fieldDT = 0,

	arrowHeight = 0,
	arrowLatDim = 30,
	arrowLonDim = 60,
	arrowEvenDistribute = true,
	arrowScale = 5,
	-- set to 0 to flatten vector field against surface
	arrowZScale = 1,
	arrowFieldNormalize = false,

	fieldLineHeight = 0,
	fieldLineLatDim = 15,
	fieldLineLonDim = 30,
	fieldLineEvenDistribute = false,
	fieldLineIter = 400,
	fieldLineIntStep = 1e-2,
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
	local chart = assert.index(allCharts, name)
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
	s = (length(B) - B3Min.x) / (B3Max.x - B3Min.x);
]],
	},
	{
		name = '|B| 2D',
		code = [[
	s = (length(B.xy) - B3Min.y) / (B3Max.y- B3Min.y);
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
	s = (B.x - BMin.x) / (BMax.x - BMin.x);
]],
	},
	{
		name = 'By (east)',
		code = [[
	s = (B.y - BMin.y) / (BMax.y - BMin.y);
]],
	},
	{
		name = 'Bz (inward)',
		code = [[
	s = (B.z - BMin.z) / (BMax.z - BMin.z);
]],
	},
	{
		name = 'div B',
		code = [[
	s = (B2.x - B2Min.x) / (B2Max.x - B2Min.x);
]],
	},
	{
		name = 'div2D B',
		code = [[
	s = (B2.y - B2Min.y) / (B2Max.y - B2Min.y);
]],
	},
	{
		name = 'curl B z',
		code = [[
	s = (B2.z - B2Min.z) / (B2Max.z - B2Min.z);
]],
	},
	{
		name = 'curl B Mag',
		code = [[
	s = (B2.w - B2Min.w) / (B2Max.w - B2Min.w);
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
				texs = {earthtex, gradTex},
			}
		end

		self.sceneobj.program = shader
		self.sceneobj.texs[2] = gradTex
		self.sceneobj.uniforms.mvProjMat = app.view.mvProjMat.ptr
		self.sceneobj.uniforms.dt = guivars.fieldDT
		self.sceneobj.uniforms.alpha = guivars.drawAlpha
		self.sceneobj.uniforms.gradScale = guivars.gradScale
		self.sceneobj.uniforms.BMin = app.BMin.s
		self.sceneobj.uniforms.BMax = app.BMax.s
		self.sceneobj.uniforms.B2Min = app.B2Min.s
		self.sceneobj.uniforms.B2Max = app.B2Max.s
		self.sceneobj.uniforms.B3Min = app.B3Min.s
		self.sceneobj.uniforms.B3Max = app.B3Max.s
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

	self.calcBCode = template(assert(path'calc_b.shader':read()), {
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

	-- lat/lon dim for surface statistics calculations
	-- used for determining ranges, tho those ranges are invalide for any other times and altitudes
	-- and in js-emulation it runs very slow
	-- so i might try to get rid of this ...
	local londim = 128	-- 1440	 -- dimension in 2D x dir / spherical phi / globe lambda dir
	local latdim = 64	-- 720	 -- dimension in 2D y dir / spherical theta / globe phi dir

	-- fbo is used for stats calcs and for field line integration
	self.fbo = GLFBO()
		:unbind()
glreport'here'

	self.quadGeomVertexShader = GLProgram.VertexShader{
		version = 'latest',
		precision = 'best',
		code = [[
layout(location=0) in vec2 vertex;
out vec2 texcoordv;
void main() {
	texcoordv = vertex;
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
]],
	}

	local piDef = 'const float M_PI = '..('%.49f'):format(math.pi)..';'

	-- used for the B2Tex, which itself is used for magnitude ranges
	-- also used for overlay obj's fragment code, for the same thing
	-- module depends M_PI, uniform vec2 latLonDim
	self.calcB2Code = [[
vec4 calcB2(vec3 plh) {
	float latdim = latLonDim.x;
	float londim = latLonDim.y;
	vec3 dphi = vec3(M_PI / londim, 0., 0.);
	vec3 dlambda = vec3(0., 2. * M_PI / latdim, 0.);
	vec3 dheight = vec3(0., 0., 1000.);

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

	return vec4(
		div_B,
		div2D_B,
		curl_B.z,
		length(curl_B)
	);
}
]]
	-- lat/lon dim for surface statistics calculations
	-- used for determining ranges, tho those ranges are invalide for any other times and altitudes
	-- TODO recalculate these when 'dt' changes
	local londim = 1440
	local latdim = 720
	local function makeLatLonFloatTex()
		return GLTex2D{
			internalFormat = gl.GL_RGBA32F,
			width = londim,
			height = latdim,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			minFilter = gl.GL_NEAREST,
			magFilter = gl.GL_NEAREST,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}:unbind()
	end

	-- fbo is used for stats calcs and for field line integration
	self.fbo = GLFBO()
		:unbind()
glreport'here'

	local piDef = '#define M_PI '..('%.49f'):format(math.pi)

	self.BTex = makeLatLonFloatTex()
	self.B2Tex = makeLatLonFloatTex()
	self.B3Tex = makeLatLonFloatTex()

	self.calcBTexSceneObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			header = piDef,
			shaders = {self.quadGeomVertexShader},
			fragmentCode = [[
uniform float dt;
]]..self.calcBCode..[[
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

glreport'here'

	-- module depends M_PI, uniform vec2 latLonDim
	self.calcB2Code = [[
vec4 calcB2(vec3 plh) {
	float latdim = latLonDim.x;
	float londim = latLonDim.y;
	vec3 dphi = vec3(M_PI / londim, 0., 0.);
	vec3 dlambda = vec3(0., 2. * M_PI / latdim, 0.);
	vec3 dheight = vec3(0., 0., 1000.);

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

	return vec4(
		div_B,
		div2D_B,
		curl_B.z,
		length(curl_B)
	);
}
]]

	self.calcB2TexSceneObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			header = piDef,
			shaders = {self.quadGeomVertexShader},
			fragmentCode = [[
uniform float dt;

]]..self.calcBCode..[[

// used in (d/dphi, d/dlambda) finite-difference resolution
uniform vec2 latLonDim;	//(latdim, londim) is (height, width)

]]..self.calcB2Code..[[

in vec2 texcoordv;
out vec4 fragColor;

void main() {
	float phi = (texcoordv.y - .5) * M_PI;			//[-pi/2, pi/2]
	float lambda = (texcoordv.x - .5) * 2. * M_PI;	//[-pi, pi]
	fragColor = calcB2(vec3(phi, lambda, 0.));
}
]],
			uniforms = {
				dt = 0,
				latLonDim = {latdim, londim},
			},
		},
		geometry = self.quadGeom,
	}

	self.calcB3TexSceneObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			shaders = {self.quadGeomVertexShader},
			fragmentCode = [[
in vec2 texcoordv;
uniform sampler2D BTex;
out vec4 fragColor;
void main() {
	vec4 B = texture(BTex, texcoordv);
	fragColor = vec4(length(B.xyz), length(B.xy), 0., 1.);
}
]],
			uniforms = {
				BTex = 0,
			},
		},
		texs = {self.BTex},
		geometry = self.quadGeom,
	}

	local GLReduce = require 'reduce'
	local reducePP = GLReduce:makePingPong{tex=self.BTex, fbo=self.fbo}
	self.minReduce = GLReduce{
		fbo = self.fbo,
		pingpong = reducePP,
		geometry = self.quadGeom,
		vertexShader = self.quadGeomVertexShader,
		gpuop = function(a,b) return 'min('..a..', '..b..')' end,
		cpuop = function(a,b,c)
			for i=0,3 do
				c.s[i] = math.min(a.s[i], b.s[i])
			end
		end
	}
	self.maxReduce = GLReduce{
		fbo = self.fbo,
		pingpong = reducePP,
		geometry = self.quadGeom,
		vertexShader = self.quadGeomVertexShader,
		gpuop = function(a,b) return 'max('..a..', '..b..')' end,
		cpuop = function(a,b,c)
			for i=0,3 do
				c.s[i] = math.max(a.s[i], b.s[i])
			end
		end
	}

	self:recalcBStats()

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
			fragmentCode =
piDef..[[

uniform float dt;

]]..self.calcBCode..[[

// used in (d/dphi, d/dlambda) finite-difference resolution
uniform vec2 latLonDim;	//(latdim, londim) is (height, width)

]]..self.calcB2Code..[[

uniform vec4 BMin, BMax;		// x y z
uniform vec4 B2Min, B2Max;		// div3d div2d curl2d curl3d
uniform vec4 B3Min, B3Max;		// mag3d mag2d

in vec2 texcoordv;
out vec4 fragColor;

uniform sampler2D earthTex;
uniform sampler2D gradTex;
uniform float alpha;
uniform float gradScale;

void main() {
	float s = .5;
	float hsvBlend = .5;

	vec3 plh = vec3(
		(texcoordv.y - .5) * M_PI,			//phi
		(texcoordv.x - .5) * 2. * M_PI,		//lambda
		0.
	);
	vec3 B = calcB(plh);
	vec4 B2 = calcB2(plh);

]]..overlay.code..[[

	fragColor = mix(
		texture(earthTex, vec2(texcoordv.x, 1. - texcoordv.y)),
		texture(gradTex, vec2(s * gradScale, .5)),
		hsvBlend);
	fragColor.a = alpha;
}
]],
			uniforms = table({
				earthTex = 0,
				gradTex = 1,
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
]]..self.calcBCode
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

	-- integrate the field lines and store as cols in a texture
	self.integrateFieldLinesSceneObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			shaders = {self.quadGeomVertexShader},
			fragmentCode = chartCode
..[[
uniform float dt;
]]..self.calcBCode..[[

// ported from WMM2020 GeomagnetismLibrary.c
// expects xyz in cartesian units earth-semimajor-axis
// output is in (radians, radians, km)
// TODO just use charts.WGS84:chartInv(x,y,z) ?  or use this there?
vec3 cartesianToLatLonWGS84(vec3 pos) {
	// lowercase wgs84_a is in km, uppercase WGS84_a is in m ... i know i need to fix this ...
	pos *= wgs84_a;	// convert from semimajor-axis units to km
	float modified_b = pos.z < 0. ? -wgs84_b : wgs84_b;
	float r = length(pos.xy);

	float e = (modified_b * pos.z - (wgs84_a * wgs84_a - modified_b * modified_b)) / (wgs84_a * r);
	float f = (modified_b * pos.z + (wgs84_a * wgs84_a - modified_b * modified_b)) / (wgs84_a * r);
	float p = (4. / 3.) * (e * f + 1.);
	float q = 2. * (e * e - f * f);
	float d = p * p * p + q * q;

	float v;
	if  (d >= 0.) {
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
		zlong += 2. * M_PI;
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

#if 0	//... and recompute the new coord's B and store the new Bmag ... wasted calcs ... I could avoid by storing both phi,lambda,hegiht and storing Bvec as we integrate,  but that's one extra texture ...
	B = calcB(phiLambdaHeightBMag.xyz);
#endif // ... or just use the previous iterations' B field magnitude, so the lines will be one vertex off in their coloring, but we'll call 1/2 the calcB()'s
	phiLambdaHeightBMag.w = length(B);
	// TODO maybe, track two tex states: pos and B, and then you can have correct indexed BMag without double the calcs

	fragColor = phiLambdaHeightBMag;
}
]],
			uniforms = {
				fieldLinePosTex = 0,
				--dt = 0,
				--dparam = guivars.fieldLineIntStep,
				--texHeight = #startCoords,
			},
		},
		geometry = self.quadGeom,
		--texs = {self.fieldLinePosTex},	-- but don't set it here, do it in the render loop
	}

	-- draw field lines
	self.fieldLineShader = GLProgram{
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
layout(location=0) in vec4 vertex;

out float BFrac;

uniform mat4 mvProjMat;
void main() {
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
]],
		{
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

	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)

	self:updateFieldLines()
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

-- if fieldDT, fieldLineLonDim, fieldLineLatDim, fieldLineEvenDistribute, fieldLineHeight change
-- then rebuild the initial spots, resize buffers, and integrate ...
function App:updateFieldLines()
	local londim = guivars.fieldLineLonDim	-- field line lon dim
	local latdim = guivars.fieldLineLatDim	-- field line lat dim
	local evenDistribute = guivars.fieldLineEvenDistribute	-- field line even distribute
	local height0 = guivars.fieldLineHeight	-- field line initial height

--[[
how to build the field lines on the GPU ...
1) pick our seed locations / fill them in on the top row of a float tex
2) run a 1 x n shader / fbo on each successive rows, integrate as we go ...
3) copy the fbo tex to an array buffer, and use it as a bunch of line strip vertexes
--]]

	self.startCoords = vector'vec4f_t'
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
				self.startCoords:emplace_back()[0]:set(phi, lambda, height0, Bmag)
			end
		end
	end
print('# field line start coords', #self.startCoords)

	-- perform integration on the GPU ...

	-- store current state here
	-- I couldve used a pingpong but why waste the memory
	self.fieldLinePosTex = GLTex2D{
		internalFormat = gl.GL_RGBA32F,	-- phi, lambda, height, |B|
		width = 1,			-- 1D along height for blitting's sake
		height = #self.startCoords,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_NEAREST,
		wrap = {
			s = gl.GL_REPEAT,
			t = gl.GL_REPEAT,
		},
		data = self.startCoords.v,
	}:unbind()

	-- use a fbo to integrate / generate our field line positions
	-- seed on the y axis so that successive lines form along the x-axis i.e. along contiguous memory to be later passed off to the draw calls
	-- (use a PBO?)
	self.fieldLineVtxsTex = GLTex2D{
		internalFormat = gl.GL_RGBA32F,	-- phi, lambda, height, |B|
		width = guivars.fieldLineIter,
		height = #self.startCoords,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_NEAREST,
		wrap = {
			s = gl.GL_REPEAT,
			t = gl.GL_REPEAT,
		},
	}:subimage{		-- upload initial state as a vertical strip
		width = 1,
		height = #self.startCoords,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		data = self.startCoords.v,
	}:unbind()

	-- copy from tex to array buffer with PBO's and glReadPixels from the FBO attached to the Tex to copy (smh...)
	self.fieldLineVertexBuf = GLArrayBuffer{
		size = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height * 4 * 4,	-- sizeof(float) * 4 components
		type = gl.GL_FLOAT,
		count = self.fieldLineVtxsTex.width * self.fieldLineVtxsTex.height,
		dim = 4,
		mode = gl.GL_STREAM_DRAW,--gl.GL_STATIC_DRAW,
	}:unbind()

	self.fieldLineSceneObj = GLSceneObject{
		program = self.fieldLineShader,
		vertexes = self.fieldLineVertexBuf,
		geometries = range(0,#self.startCoords-1):mapi(function(i)
			return {
				mode = gl.GL_LINE_STRIP,
				count = guivars.fieldLineIter,
				offset = i * guivars.fieldLineIter,
			}
		end),
	}

	self:integrateFieldLines()
end

-- if fieldDT changes then reintegrate
-- TODO fieldLineHeight ... but that'd mean repopulating the initial z channel
function App:integrateFieldLines()
	-- do the integration using FBO's onto our texture

	self.fbo:bind()
		:setColorAttachmentTex2D(self.fieldLineVtxsTex.id)
	local res, err = self.fbo.check()
	if not res then print(err) end
	gl.glReadBuffer(gl.GL_COLOR_ATTACHMENT0)

	self.integrateFieldLinesSceneObj.program:use()
	self.integrateFieldLinesSceneObj.program:setUniforms{
		dt = guivars.fieldDT,
		dparam = guivars.fieldLineIntStep,
		texHeight = self.fieldLineVtxsTex.height,	-- #startCoords
	}

	-- now iteratively draw to FBO next strip over, reading from the current state strip as we go
	self.fieldLinePosTex:bind()
	for i=1,guivars.fieldLineIter-1 do
		-- draw from fieldLinePosTex into fieldLineVtxsTex
		gl.glViewport(i, 0, 1, self.fieldLineVtxsTex.height)
		self.integrateFieldLinesSceneObj:draw()
		-- copy from fbo fieldLineVtxsTex back into fieldLinePosTex for the next iteration
		gl.glCopyTexSubImage2D(self.fieldLinePosTex.target, 0, 0, 0, i, 0, 1, self.fieldLineVtxsTex.height)
	end
	-- final pass, copy the initial state back into fieldLinePosTex, so the next integration pass will use it as the initial state
	gl.glCopyTexSubImage2D(self.fieldLinePosTex.target, 0, 0, 0, 0, 0, 1, self.fieldLineVtxsTex.height)

	self.fieldLinePosTex:unbind()

	-- now use a PBO to copy the tex into a vertex attribute

	self.fieldLineVertexBuf:bind(gl.GL_PIXEL_PACK_BUFFER)
	gl.glViewport(0, 0, self.fieldLineVtxsTex.width, self.fieldLineVtxsTex.height)
	-- read from bound FBO (attached to tex) to PBO (attached to buffer)
	gl.glReadPixels(0, 0, self.fieldLineVtxsTex.width, self.fieldLineVtxsTex.height, gl.GL_RGBA, gl.GL_FLOAT, nil)
	self.fieldLineVertexBuf:unbind(gl.GL_PIXEL_PACK_BUFFER)

	self.fbo:unbind()

	gl.glReadBuffer(gl.GL_BACK)
glreport'here'
end

-- this goes slow right now
-- doing a min/max reduce on my hydro opencl project is muuuuuch faster
-- in fact, am I doing reduce on my webgl cfd project?
-- hmm, how to improve and still be glsl compat
-- for one, I could just reduce the variable being drawn in overlay, and not reduce all of them
function App:recalcBStats()
	self.calcBTexSceneObj.uniforms.dt = guivars.fieldDT
	self.fbo:draw{
		viewport = {0, 0, self.BTex.width, self.BTex.height},
		dest = self.BTex,
		callback = function()
			self.calcBTexSceneObj:draw()
		end,
	}

	self.calcB2TexSceneObj.uniforms.dt = guivars.fieldDT
	self.fbo:draw{
		viewport = {0, 0, self.BTex.width, self.BTex.height},
		dest = self.B2Tex,
		callback = function()
			self.calcB2TexSceneObj:draw()
		end,
	}

	self.calcB3TexSceneObj.uniforms.dt = guivars.fieldDT
	self.fbo:draw{
		viewport = {0, 0, self.BTex.width, self.BTex.height},
		dest = self.B3Tex,
		callback = function()
			self.calcB3TexSceneObj:draw()
		end,
	}

	timer('generating stats', function()
-- [[
		self.BMin = self.minReduce(self.BTex)
		self.BMax = self.maxReduce(self.BTex)
--print('reduce B min', BMin, 'max', BMax)
		self.B2Min = self.minReduce(self.B2Tex)
		self.B2Max = self.maxReduce(self.B2Tex)
--print('reduce B2 min', B2Min, 'max', B2Max)
		self.B3Min = self.minReduce(self.B3Tex)
		self.B3Max = self.maxReduce(self.B3Tex)
--print('reduce B3 min', B3Min, 'max', B3Max)
--]]

--[[ stats should look like for wmm2020 dt=0
x = {min = -16735.287109375, max = 41797.078125, avg = 17984.161021002, sqavg = 460231098.77605, stddev = 11696.198149258, count = 1036800},
y = {min = -17571.607421875, max = 16772.9140625, avg = 0.00060528925769373, sqavg = 42031905.017686, stddev = 6483.2017566698, count = 1036800},
z = {min = -66933.3984375, max = 60970.375, avg = 1180.7201683604, sqavg = 1732190539.499, stddev = 41602.841722447, count = 1036800},
div = {min = -4.3255414962769, max = 4.7061634063721, avg = 0.0019280310092957, sqavg = 0.044430527339689, stddev = 0.2107766828568, count = 1036800},
div2d = {min = -4.2783694267273, max = 4.6579914093018, avg = 0.0029287106643504, sqavg = 0.039302149729322, stddev = 0.19822606383411, count = 1036800},
curlZ = {min = -104.66393280029, max = 103.64836120605, avg = 1.1786829213153e-08, sqavg = 7.8592639803675, stddev = 2.8034378859478, count = 1036800},
curlMag = {min = 0.00051679083844647, max = 104.75435638428, avg = 0.16494477450588, sqavg = 7.9169097917121, stddev = 2.8088615154677, count = 1036800},
mag = {min = 22232.017209054, max = 66990.328957622, avg = 45853.59896298, sqavg = 2234453543.2928, stddev = 11484.816299572, count = 1036800},
mag2d = {min = 29.5594549176, max = 41800.698442297, avg = 20242.527057979, sqavg = 502263003.79371, stddev = 9617.85330002, count = 1036800},
--]]
	end)
--]=]
glreport'here'
end

--[=[
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

	gl.glEnable(gl.GL_DEPTH_TEST)

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
		self.arrowFieldSceneObj.uniforms.BMagMax = self.B3Max.x
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
		self.fieldLineSceneObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
		self.fieldLineSceneObj.uniforms.BMag = {self.B3Min.x, self.B3Max.x}

		-- TODO if we're using weights then why do we need the 'chartIs3D' flag?
		self.fieldLineSceneObj.uniforms.weight_Equirectangular = guivars.geomIndex == chartIndexForName.Equirectangular and 1 or 0
		self.fieldLineSceneObj.uniforms.weight_Azimuthal_equidistant = guivars.geomIndex == chartIndexForName['Azimuthal equidistant'] and 1 or 0
		self.fieldLineSceneObj.uniforms.weight_Mollweide = guivars.geomIndex == chartIndexForName.Mollweide and 1 or 0
		self.fieldLineSceneObj.uniforms.weight_WGS84 = guivars.geomIndex == chartIndexForName.WGS84 and 1 or 0
		self.fieldLineSceneObj.uniforms.chartIs3D = chart.is3D or false

		self.fieldLineSceneObj:draw()
	end

	gl.glEnable(gl.GL_BLEND)

	chart:draw(self, shader, gradTex)

	gl.glDisable(gl.GL_BLEND)
	gl.glDisable(gl.GL_DEPTH_TEST)

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

	if ig.luatableInputFloat('field line init altitude', guivars, 'fieldLineHeight')	-- techniically this should just trigger a re-integration, not a realloc
	or ig.luatableInputInt('field line lat dim', guivars, 'fieldLineLatDim')
	or ig.luatableInputInt('field line lon dim', guivars, 'fieldLineLonDim')
	or ig.luatableCheckbox('field lines even distribute', guivars, 'fieldLineEvenDistribute')
	or ig.luatableInputInt('field line iter', guivars, 'fieldLineIter')
	then
		self:updateFieldLines()
	end

	if ig.luatableInputFloat('field line int step', guivars, 'fieldLineIntStep') then
		self:integrateFieldLines()
	end

	-- how linear are the g and h coeffs?
	-- can I just factor out the dt?
	--ig.luatableInputFloat('time from '..wmm.epoch, guivars, 'fieldDT')
	if ig.luatableSliderFloat('years from '..tostring(wmm.epoch), guivars, 'fieldDT', -50, 50) then
		self.lastUpdateDTTime = timer.getTime()
		self:integrateFieldLines()
	end

	-- only once the dt var has finished updating ...
	if self.lastUpdateDTTime
	and timer.getTime() - self.lastUpdateDTTime > 1
	then
		self.lastUpdateDTTime = nil
		self:recalcBStats()
	end
end

return App():run()

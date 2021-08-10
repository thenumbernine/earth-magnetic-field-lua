#!/usr/bin/env luajit
local ffi = require 'ffi'
local vec3f = require 'vec-ffi.vec3f'
local template = require 'template'
local gl = require 'gl'
local glreport = require 'gl.report'
local ig = require 'ffi.imgui'
local GLTex2D = require 'gl.tex2d'
local GLProgram = require 'gl.program'
local glCallOrRun = require 'gl.call'
require 'ext'


-- Sets WGS-84 parameters
local wgs84 = {}
wgs84.a = 6378.137 --semi-major axis of the ellipsoid in 
wgs84.b = 6356.7523142 --semi-minor axis of the ellipsoid in 
wgs84.fla = 1 / 298.257223563 -- flattening 
wgs84.eps = math.sqrt(1 - (wgs84.b * wgs84.b) / (wgs84.a * wgs84.a)) --first eccentricity 
wgs84.epssq = wgs84.eps * wgs84.eps --first eccentricity squared 
wgs84.re = 6371.2 -- Earth's radius 



local londim = 1440	-- dimension in 2D x dir / spherical phi / globe lambda dir
local latdim = 720	-- dimension in 2D y dir / spherical theta / globe phi dir
local HeightAboveEllipsoid = 0		-- compute at z=0 for now
local year = 2020

-- compute a 2D grid of the field
local Bdata = ffi.new('vec3f_t[?]', londim * latdim)
-- TODO later -- compute a 3D grid

local BStat = require 'stat.set'('mag', 'x', 'y', 'z')

-- cache numbers
local fn = ('bfield_year=%d_londim=%d_latdim=%d.f32'):format(year, londim, latdim)
if os.fileexists(fn) then
	local data = file[fn]
	local s = ffi.cast('char*', data)
	local f = ffi.cast('vec3f_t*', s)
	assert(#data == londim * latdim * ffi.sizeof'vec3f_t')

	ffi.copy(ffi.cast('char*', Bdata), f, #data)

	for j=0,latdim-1 do
		for i=0,londim-1 do
			local e = i + londim * j
			local B = Bdata[e]
			local Bmag = B:length()
			BStat:accum(B:length(), B.x, B.y, B.z)
		end
	end
else
	local lines = file['wmm.cof']:trim():split'\n'
	local header = lines:remove(1):trim():split'%s+'
	assert(#header == 3)
	print('model epoch', header[1])
	print('model name', header[2])
	print('date of release', header[3])
	assert(lines:remove() == '999999999999999999999999999999999999999999999999')
	assert(lines:remove() == '999999999999999999999999999999999999999999999999')

	local wmm = {}
	for _,line in ipairs(lines) do
		local parts = line:trim():split'%s+'
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
			gt = parts[5],
			ht = parts[6],
		}
		-- n runs 1-12
		-- m runs 0-n
	end

	local nMax = #wmm


	for j=0,latdim-1 do
		local phi = ((j+.5)/latdim * 2 - 1) * 90	-- spherical theta
		phi = math.rad(phi)
			
		-- begin MAG_GeodeticToSpherical
		local CosLat = math.cos(phi)
		local SinLat = math.sin(phi)

		-- convert from geodetic WGS-84 to spherical coordiantes
		local rc = wgs84.a / math.sqrt(1 - wgs84.epssq * SinLat * SinLat)
		local xp = (rc + HeightAboveEllipsoid) * CosLat
		local zp = (rc * (1 - wgs84.epssq) + HeightAboveEllipsoid) * SinLat
		
		-- spherical results:
		local r = math.sqrt(xp * xp + zp * zp)
		local phig = math.asin(zp / r) -- geocentric latitude 
		-- longitude is the same 
		-- end MAG_GeodeticToSpherical

		for i=0,londim-1 do
			local lambda = ((i+.5)/londim * 2 - 1) * 180	-- spherical phi
			lambda = math.rad(lambda)
	   
			-- begin MAG_Geomag
			-- begin MAG_ComputeSphericalHarmonicVariables

			local cos_lambda = math.cos(lambda)
			local sin_lambda = math.sin(lambda)
			
			local RelativeRadiusPower = {}	-- 0-nmake
			RelativeRadiusPower[0] = (wgs84.re / r) * (wgs84.re / r)
			for n=1,nMax do
				RelativeRadiusPower[n] = RelativeRadiusPower[n-1] * (wgs84.re / r)
			end

			local cos_mlambda = {[0] = 1, cos_lambda}
			local sin_mlambda = {[0] = 0, sin_lambda}

			-- looks like exp(i lambda)
			for m=2,nMax do
				cos_mlambda[m] = cos_mlambda[m-1] * cos_lambda - sin_mlambda[m-1] * sin_lambda
				sin_mlambda[m] = cos_mlambda[m-1] * sin_lambda + sin_mlambda[m-1] * cos_lambda
			end
			
			-- end MAG_ComputeSphericalHarmonicVariables
			-- begin MAG_AssociatedLegendreFunction

			local sin_phi = math.sin(phig)	-- why convert to degrees?
		
			-- begin MAG_PcupLow
			
			local x = sin_phi
			local Pcup = {[0] = 1}
			local dPcup = {[0] = 0}

			-- sin (geocentric latitude) - sin_phi
			local z = math.sqrt((1 - x) * (1 + x));

			local NumTerms = ((nMax + 1) * (nMax + 2) / 2)
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
							k = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3))
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
							wmm_n_m.g * cos_mlambda[m] 
							+ wmm_n_m.h * sin_mlambda[m]
						)
						* (n + 1) * Pcup[index]

					--		  1 nMax  (n+2)    n     m            m           m
					--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
					--				   n=1             m=0   n            n           n  */
					-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
					By = By + (
						RelativeRadiusPower[n]
						* (
							wmm_n_m.g * sin_mlambda[m] 
							- wmm_n_m.h * cos_mlambda[m]
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
							wmm_n_m.g * cos_mlambda[m]
							+ wmm_n_m.h * sin_mlambda[m]
						)
						* dPcup[index]
				end
			end

			local cos_phi = math.cos(phig)
			if math.abs(cos_phi) > 1e-10 then
				By = By / cos_phi;
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
					local schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt( (n * 2) /  (n + 1))
					local schmidtQuasiNorm1 = schmidtQuasiNorm2
					if n == 1 then
						PcupS[n] = PcupS[n-1]
					else
						local k =  (((n - 1) * (n - 1)) - 1) /  ((2 * n - 1) * (2 * n - 3))
						PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2]
					end

					--		  1 nMax  (n+2)    n     m            m           m
					--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
					--				   n=1             m=0   n            n           n  */
					-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
					By = By + 
						RelativeRadiusPower[n] *
						(
							wmm_n_m.g * sin_mlambda[1] -
							wmm_n_m.h * cos_mlambda[1])
							* PcupS[n] * schmidtQuasiNorm3
				end

				-- end MAG_SummationSpecial	
			end
			
			-- end MAG_Summation 
			-- end MAG_Geomag

			local e = i + londim * j
			
			Bdata[e]:set(Bx, By, Bz)
			local B = Bdata[e]
			BStat:accum(B:length(), B.x, B.y, B.z)
		end
	end

	file[fn] = ffi.string(ffi.cast('char*', Bdata), londim * latdim * ffi.sizeof'vec3f_t')
end

print('BStat')
print(BStat)

local earthtex
local Btex

local App = class(require 'glapp.orbit'(require 'imguiapp'))

App.title = 'EM field' 

local geomIndex = ffi.new('int[1]', 0)
local overlayIndex = ffi.new('int[1]', 1)
local gradientIndex = ffi.new('int[1]', 0)

local geoms = {
	{
		name = '2D',
		draw = function()
			gl.glBegin(gl.GL_QUADS)
			gl.glTexCoord2f(0, 0)	gl.glVertex2f(-2, -1)
			gl.glTexCoord2f(1, 0)	gl.glVertex2f(2, -1)
			gl.glTexCoord2f(1, 1)	gl.glVertex2f(2, 1)
			gl.glTexCoord2f(0, 1)	gl.glVertex2f(-2, 1)
			gl.glEnd()
		end,
	},
	{
		name = '3D',
		draw = function(self)
			self.list = self.list or {}
			glCallOrRun(self.list, function()
				local HeightAboveEllipsoid = 0
				local jres = 120
				local ires = 60
				for ibase=0,ires-1 do
					gl.glBegin(gl.GL_TRIANGLE_STRIP)
					for j=0,jres do
						local v = j/jres
						local lambda = math.rad((v * 2 - 1) * 180)
						local cos_lambda = math.cos(lambda)
						local sin_lambda = math.sin(lambda)
						for iofs=0,1 do
							local i = ibase + iofs
							local u = i/ires
							local phi = math.rad((u * 2 - 1) * 90)
							local CosLat = math.cos(phi)
							local SinLat = math.sin(phi)
							
							local rc = wgs84.a / math.sqrt(1 - wgs84.epssq * SinLat * SinLat)
							local xp = (rc + HeightAboveEllipsoid) * CosLat
							local zp = (rc * (1 - wgs84.epssq) + HeightAboveEllipsoid) * SinLat
							local r = math.sqrt(xp * xp + zp * zp)
							local phig = math.asin(zp / r)
							-- TODO check domain and see if you can just use zp/r and sqrt(1-sin_phi^2)
							local cos_phi = math.cos(phig)
							local sin_phi = math.sin(phig)
						
							r = r / wgs84.a
							local x = r * math.cos(phig) * cos_lambda
							local y = r * math.cos(phig) * sin_lambda
							local z = r * math.sin(phig)

							gl.glTexCoord2f(v, u)
							gl.glVertex3f(x, y, z)
						end
					end
					gl.glEnd()
				end
			end)
		end,
	},
}

local gradients = {
	{
		name = 'rainbow',
		gen = function()
			return require 'gl.hsvtex'(256)
		end,
	},
	{
		name = 'B&W',
		gen = function()
			local n = 256
			-- why is this intermittantly misaligned?
			local image = require 'image'(
				n, 1, 4, 'unsigned char', function(u)
					if bit.band(u, 15) == 0 then
						return 255,255,255,255
					else
						return 0,0,0,0
					end
				end
			)
			return require 'gl.tex1d'{
				--image = image,
				data = image.buffer,
				width = image.width,

				internalFormat = gl.GL_RGBA,
				format = gl.GL_RGBA,
				type = gl.GL_UNSIGNED_BYTE,
				minFilter = gl.GL_LINEAR,
				--magFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
				magFilter = gl.GL_LINEAR,
				generateMipmap = true,
			}
		end,
	},
}

local overlays = {
	{
		name = 'Earth',
		code = [[
	alpha = 0.;
]],
	},
	{
		name = '|B|', 
		code = [[
	s = (length(B) - BMagMin) / (BMagMax - BMagMin);
]],
	},
	{
		name = 'arg B',
		code = [[
	s = atan(B.y, B.x) / (2. * M_PI) + .5;
]],
	},
	{
		name = 'Bx',
		code = [[
	s = (B.x - BxMin) / (BxMax - BxMin);
]],
	},
	{
		name = 'By',
		code = [[
	s = (B.y - ByMin) / (ByMax - ByMin);
]],
	},
	{
		name = 'Bz',
		code = [[
	s = (B.z - BzMin) / (BzMax - BzMin);
]],
	},
}



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

	Btex = GLTex2D{
		internalFormat = gl.GL_RGBA32F,
		width = londim,
		height = latdim,
		format = gl.GL_RGB,
		type = gl.GL_FLOAT,
		data = ffi.cast('char*', Bdata),
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_NEAREST,
	}

	for _,overlay in ipairs(overlays) do
		overlay.shader = GLProgram{
			vertexCode = [[
varying vec2 tc;
void main() {
	tc = gl_MultiTexCoord0.xy;
	gl_Position = ftransform();
}
]],
			fragmentCode = template([[
#define BMagMin <?=clnumber(BStat.mag.min)?>
#define BMagMax <?=clnumber(BStat.mag.max)?>
#define BxMin <?=clnumber(BStat.x.min)?>
#define BxMax <?=clnumber(BStat.x.max)?>
#define ByMin <?=clnumber(BStat.y.min)?>
#define ByMax <?=clnumber(BStat.y.max)?>
#define BzMin <?=clnumber(BStat.z.min)?>
#define BzMax <?=clnumber(BStat.z.max)?>
#define M_PI <?=('%.49f'):format(math.pi)?>

varying vec2 tc;

uniform sampler2D earthtex;
uniform sampler2D Btex;
uniform sampler1D hsvtex;

void main() {
	float s = .5;
	float alpha = .5;
	vec3 B = texture2D(Btex, tc).rgb;
	<?=overlay.code?>
	gl_FragColor = mix(
		texture2D(earthtex, vec2(tc.x, 1. - tc.y)),
		texture1D(hsvtex, s),
		alpha);
}
]], 		{
				overlay = overlay,
				BStat = BStat,
				clnumber = require 'cl.obj.number',
			}),
			uniforms = {
				earthtex = 0,
				Btex = 1,
				hsvtex = 2,
			},
		}
		overlay.shader:useNone() 
	end

	if self.view then
		self.view.ortho = true
		self.view.orthoSize = 2
	end
	gl.glClearColor(0,0,0,0)

	gl.glEnable(gl.GL_DEPTH_TEST)
end

function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local shader = assert(overlays[tonumber(overlayIndex[0])+1]).shader
	local gradtex = assert(gradients[tonumber(gradientIndex[0])+1]).tex
	local geom = assert(geoms[tonumber(geomIndex[0])+1])

	shader:use()
	
	earthtex:bind(0)
	Btex:bind(1)
	gradtex:bind(2)

	geom:draw()
	
	gradtex:unbind(2)
	Btex:unbind(1)
	earthtex:unbind(0)
	
	shader:useNone()
	
	glreport'here'

	--render gui
	App.super.update(self, ...)
end

local bool = ffi.new('bool[1]', false)

function App:updateGUI()
	ig.igText'geom'
	for i,geom in ipairs(geoms) do
		ig.igRadioButtonIntPtr(geom.name, geomIndex, i-1)
	end

	ig.igSeparator()
	
	ig.igText'overlay'
	for i,overlay in ipairs(overlays) do
		ig.igRadioButtonIntPtr(overlay.name, overlayIndex, i-1)
	end
	
	ig.igSeparator()
	
	ig.igText'gradient'
	for i,grad in ipairs(gradients) do
		ig.igRadioButtonIntPtr(grad.name, gradientIndex, i-1)
	end

	bool[0] = self.view.ortho
	if ig.igCheckbox('ortho', bool) then
		self.view.ortho = bool[0]
	end
end

App():run()

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




-- load wmm
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



local function calcB(phi, lambda, HeightAboveEllipsoid)

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
			local schmidtQuasiNorm3 = schmidtQuasiNorm2 * math.sqrt( (n * 2) /  (n + 1))
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

	return Bx, By, Bz
end






-- compute a 2D grid of the field
local Bdata = ffi.new('vec3f_t[?]', londim * latdim)
-- TODO later -- compute a 3D grid

local BStat = require 'stat.set'('mag', 'x', 'y', 'z', 'mag2d')

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
			BStat:accum(B:length(), B.x, B.y, B.z, math.sqrt(B.x*B.x + B.y*B.y))
		end
	end
else
	for j=0,latdim-1 do
		local phi = math.rad(((j+.5)/latdim * 2 - 1) * 90)	-- spherical theta
		for i=0,londim-1 do
			local lambda = math.rad(((i+.5)/londim * 2 - 1) * 180)	-- spherical phi
			local Bx, By, Bz = calcB(phi, lambda, HeightAboveEllipsoid)
			local e = i + londim * j
			Bdata[e]:set(Bx, By, Bz)
			local B = Bdata[e]
			BStat:accum(B:length(), B.x, B.y, B.z, math.sqrt(B.x*B.x + B.y*B.y))
		end
	end

	file[fn] = ffi.string(ffi.cast('char*', Bdata), londim * latdim * ffi.sizeof'vec3f_t')
end


-- latitude = phi in [-90,90]
-- longitude = lambda in [-180,180]
local function latLonToCartesianWGS84(phi, lambda, height)
	local cosLambda = math.cos(lambda)
	local sinLambda = math.sin(lambda)

	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)
	
	local rCart = wgs84.a / math.sqrt(1 - wgs84.epssq * sinPhi * sinPhi)
	local xp = (rCart + height) * cosPhi
	local zp = (rCart * (1 - wgs84.epssq) + height) * sinPhi
	
	local r2D = math.sqrt(xp * xp + zp * zp)
	local sinPhiSph = zp / r2D
	local cosPhiSph = math.sqrt(1 - sinPhiSph * sinPhiSph)

	r2D = r2D / wgs84.a
	local x = r2D * cosPhiSph * cosLambda
	local y = r2D * cosPhiSph * sinLambda
	local z = r2D * sinPhiSph

	return x, y, z
end

--[=[
do
	local sym = require 'symmath'
	local var = sym.var
	local lambda = var'lambda'
	local phi = var'phi'
	local height = var'height'
	local wgs84_a = sym.var'wgs84.a'
	local wgs84_epssq = sym.var'wgs84.epssq'
	
	local cosLambda = sym.cos(lambda)
	local sinLambda = sym.sin(lambda)

	local cosPhi = sym.cos(phi)
	local sinPhi = sym.sin(phi)
	
	local rCart = wgs84_a / sym.sqrt(1 - wgs84_epssq * sinPhi * sinPhi)
	local xp = (rCart + height) * cosPhi
	local zp = (rCart * (1 - wgs84_epssq) + height) * sinPhi
	
	local r2D = sym.sqrt(xp * xp + zp * zp)
	local phiSph = sym.asin(zp / r2D)
	-- TODO check domain and see if you can just use zp/r2D and sqrt(1-sin_phi^2)
	local cosPhiSph = sym.cos(phiSph)
	local sinPhiSph = sym.sin(phiSph)

	r2D = r2D / wgs84_a
	local x = r2D * cosPhiSph * cosLambda
	local y = r2D * cosPhiSph * sinLambda
	local z = r2D * sinPhiSph


	local chart = sym.Array(x,y,z)

chart = chart:prune()
chart = chart:tidy()
--chart() takes forever

--print'chart:'
--print(sym.export.Lua:toFunc{func='latLonToCartesianWGS84', output={table.unpack(chart)}, input={phi, lambda, height}})

-- TODO how to evaluate only 

-- TODO divide by r cos(theta)
local Bx = sym.Array(
	x:diff(phi):prune():tidy(),
	y:diff(phi):prune():tidy(),
	z:diff(phi):prune():tidy()
)
--[[ not enough memory
BxLen = Bx:normSq():sqrt():prune():tidy()
Bx = (Bx / BxLen):prune():tidy()
--]]
-- TODO divide by r
local By = sym.Array(
	x:diff(lambda):prune():tidy(),
	y:diff(lambda):prune():tidy(),
	z:diff(lambda):prune():tidy()
)
--[[ not enough memory
ByLen = By:normSq():sqrt():prune():tidy()
By = (By / ByLen):prune():tidy()
--]]
local Bz = sym.Array(
	-x:diff(height):prune():tidy(),
	-y:diff(height):prune():tidy(),
	-z:diff(height):prune():tidy()
)
--[[ not enough memory
BzLen = Bz:normSq():sqrt():prune():tidy()
Bz = (Bz / BzLen):prune():tidy()
--]]

print'derivs:'
print(sym.export.Lua:toFunc{func='latLonToCartesianTangentSpaceWGS84', output={
	Bx[1], Bx[2], Bx[3],
	By[1], By[2], By[3],
	Bz[1], Bz[2], Bz[3],
}, input={phi, lambda, height}})

os.exit()
--	Bx = d/dphi map
--	By = d/dlambda map
--	Bz = -d/dheight map
end
--]=]

--[[
latitude = phi in [-90,90]
longitude = lambda in [-180,180]
Bx = d/dphi points north
By = d/dlambda points east
Bz = -d/dheight points inwards
--]]
local function latLonToCartesianTangentSpaceWGS84(phi, lambda, height)
	local cosLambda = math.cos(lambda)
	local sinLambda = math.sin(lambda)

	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)
	local dphi_cosPhi = -sinPhi
	local dphi_sinPhi = cosPhi

	local rCart = wgs84.a / math.sqrt(1 - wgs84.epssq * sinPhi * sinPhi)
	local dphi_rCart = wgs84.a / math.sqrt(1 - wgs84.epssq * sinPhi * sinPhi)^3 * wgs84.epssq * sinPhi * dphi_sinPhi

	local rCart_over_a = 1 / math.sqrt(1 - wgs84.epssq * sinPhi * sinPhi)

	local xp = (rCart + height) * cosPhi
	local dphi_xp = dphi_rCart * cosPhi + (rCart + height) * dphi_cosPhi
	local dheight_xp = cosPhi
	
	local xp_over_a = (rCart_over_a + height / wgs84.a) * cosPhi

	local zp = (rCart * (1 - wgs84.epssq) + height) * sinPhi
	local dphi_zp = (dphi_rCart * (1 - wgs84.epssq)) * sinPhi + (rCart * (1 - wgs84.epssq) + height) * dphi_sinPhi
	local dheight_zp = sinPhi

	local zp_over_a = (rCart_over_a * (1 - wgs84.epssq) + height / wgs84.a) * sinPhi

	local r2D = math.sqrt(xp * xp + zp * zp)
	local dphi_r2D = (xp * dphi_xp + zp * dphi_zp) / r2D
	local dheight_r2D = (xp * dheight_xp + zp * dheight_zp) / r2D
	
	local r2D_over_a = math.sqrt(xp_over_a * xp_over_a + zp_over_a * zp_over_a)
	local dphi_r2D_over_a = (xp_over_a * dphi_xp + zp_over_a * dphi_zp) / r2D

	local sinPhiSph = zp / r2D
	local dphi_sinPhiSph = (dphi_zp * r2D - zp * dphi_r2D) / (r2D * r2D)
	local dheight_sinPhiSph = (dheight_zp * r2D - zp * dheight_r2D) / (r2D * r2D)

	local cosPhiSph = math.sqrt(1 - sinPhiSph * sinPhiSph)
	--d/du sqrt(1 - x^2) = -x/sqrt(1 - x^2) dx/du
	local dphi_cosPhiSph = -sinPhi / cosPhiSph * dphi_sinPhiSph
	local dheight_cosPhiSph = -sinPhi / cosPhiSph * dheight_sinPhiSph

	--local x = r2D * cosPhiSph / wgs84.a * cosLambda
	--local y = r2D * cosPhiSph / wgs84.a * sinLambda
	--local z = r2D * sinPhiSph / wgs84.a
	
	local dphi_x = (dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * cosLambda
	local dphi_y = (dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * sinLambda
	local dphi_z = (dphi_r2D_over_a * sinPhiSph + r2D_over_a * dphi_sinPhiSph)
	
	local dlambda_x = -sinLambda
	local dlambda_y = cosLambda
	
	local dheight_x = (dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * cosLambda
	local dheight_y = (dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * sinLambda
	local dheight_z = (dheight_r2D * sinPhiSph + r2D * dheight_sinPhiSph)

	return
		vec3f(dphi_x, dphi_y, dphi_z),
		vec3f(dlambda_x, dlambda_y, 0),
		vec3f(-dheight_x, -dheight_y, -dheight_z)
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
		chart = function(phi, lambda, height) 
			return lambda*2/math.pi, phi*2/math.pi, height 
		end,
		basis = function(geom, phi, lambda, height)
			-- Bx is north, By is east, Bz is down ... smh
			return 
				vec3f(0, 1, 0),
				vec3f(1, 0, 0),
				vec3f(0, 0, -1)
		end,
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
		chart = latLonToCartesianWGS84,
		basis = function(geom, phi, lambda, height)
			return latLonToCartesianTangentSpaceWGS84(phi, lambda, height)
		end,
		draw = function(self)
			local HeightAboveEllipsoid = 0
			local jres = 120
			local ires = 60
			self.list = self.list or {}
			glCallOrRun(self.list, function()
				for ibase=0,ires-1 do
					gl.glBegin(gl.GL_TRIANGLE_STRIP)
					for j=0,jres do
						local v = j/jres
						local lambda = math.rad((v * 2 - 1) * 180)
						for iofs=0,1 do
							local i = ibase + iofs
							local u = i/ires
							local phi = math.rad((u * 2 - 1) * 90)
							
							local x,y,z = self.chart(phi, lambda, HeightAboveEllipsoid)
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
		name = 'stripes',
		gen = function()
			local n = 256
			-- why is this intermittantly misaligned?
			local image = require 'image'(
				n, 1, 4, 'unsigned char', function(u)
					local s = (u+.5)/n
					if bit.band(u, 15) == 0 then
						return 255*(1-s),127,255*s,255
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
				wrap = {
					s = gl.GL_CLAMP_TO_EDGE,
					t = gl.GL_CLAMP_TO_EDGE,
				},
				generateMipmap = true,
			}
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
		minFilter = gl.GL_LINEAR,
		magFilter = gl.GL_LINEAR,
		--magFilter = gl.GL_LINEAR_MIPMAP_LINEAR,	-- not working
		generateMipmap = true,
		wrap = {
			s = gl.GL_REPEAT,
			t = gl.GL_REPEAT,
		},
	}
	glreport'here'

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
#define BMag2DMin <?=clnumber(BStat.mag2d.min)?>
#define BMag2DMax <?=clnumber(BStat.mag2d.max)?>
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
uniform float alpha;

void main() {
	float s = .5;
	float hsvBlend = .5;
	vec3 B = texture2D(Btex, tc).rgb;
	<?=overlay.code?>
	gl_FragColor = mix(
		texture2D(earthtex, vec2(tc.x, 1. - tc.y)),
		texture1D(hsvtex, s),
		hsvBlend);
	gl_FragColor.a = alpha;
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
				alpha = 1,
			},
		}
		overlay.shader:useNone() 
	end

	if self.view then
		self.view.ortho = true
		self.view.pos.z = 2
		self.view.orthoSize = 2
	end
	gl.glClearColor(0,0,0,0)

	gl.glCullFace(gl.GL_FRONT_FACE)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
end

local function degmintofrac(deg, min, sec)
	return deg + (1/60) * (min + (1/60) * sec)
end

local function drawReading(info)
	local height = 1e-6
	local lat = degmintofrac(table.unpack(info.lat))
	local lon = degmintofrac(table.unpack(info.lon))

	gl.glColor3f(1,1,0)
	
	-- TODO geom should be transforms, and apply to all geom rendered
	gl.glPointSize(5)
	gl.glBegin(gl.GL_POINTS)

	-- reading #1
	gl.glVertex3f(
		info.geom.chart(
			math.rad(lat),
			math.rad(lon),
			height
		)
	)
	
	gl.glEnd()
	gl.glPointSize(1)

	-- draw north vs magnetic north heading
	gl.glBegin(gl.GL_LINES)
	
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

local function drawVectorField(geom)

	local arrow = {
		{-.5, 0.},
		{.5, 0.},
		{.2, .3},
		{.5, 0.},
		{.2, -.3},
		{.5, 0.},
	}

	local HeightAboveEllipsoid = 0
	local jres = 60
	local ires = 30
	local scale = 10 / (BStat.mag.max * jres)
--	geom.list = geom.list or {}
--	glCallOrRun(geom.list, function()
	gl.glBegin(gl.GL_LINES)
	for i=0,ires do
		local u = i/ires
		local phi = math.rad((u * 2 - 1) * 90)
		for j=0,jres do
			local v = j/jres
			local lambda = math.rad((v * 2 - 1) * 180)
			local x,y,z = geom.chart(phi, lambda, 0)
			gl.glTexCoord2f(u,v)
			-- TODO calc B here, and add here
			local Bx, By, Bz = calcB(phi, lambda, HeightAboveEllipsoid)
			Bx = Bx * scale	-- north component
			By = By * scale	-- east component
			Bz = Bz * scale	-- down component
			
			local BxBasis, ByBasis, BzBasis = geom:basis(phi, lambda, HeightAboveEllipsoid)

			--[[ verify the basis is correct
			scale = 3 / jres
			gl.glVertex3f(x, y, z)	gl.glVertex3f(x + BxBasis.x * scale, y + BxBasis.y * scale, z + BxBasis.z * scale)
			gl.glVertex3f(x, y, z)	gl.glVertex3f(x + ByBasis.x * scale, y + ByBasis.y * scale, z + ByBasis.z * scale)
			gl.glVertex3f(x, y, z)	gl.glVertex3f(x + BzBasis.x * scale, y + BzBasis.y * scale, z + BzBasis.z * scale)
			--]]
			-- [[ draw our arrow
			local B = BxBasis * Bx + ByBasis * By + BzBasis * Bz

			for _,q in ipairs(arrow) do
				local s, t = q[1], q[2]
				gl.glVertex3f(
					x + B.x * s - B.y * t,
					y + B.y * s + B.x * t,
					z + B.z * s)
			end
			--]]
		end
	end
	gl.glEnd()
--	end)
end


local drawAlpha = ffi.new('float[1]', 1)
local doDrawVectorField = false

function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local shader = assert(overlays[tonumber(overlayIndex[0])+1]).shader
	local gradtex = assert(gradients[tonumber(gradientIndex[0])+1]).tex
	local geom = assert(geoms[tonumber(geomIndex[0])+1])

	shader:use()

	earthtex:bind(0)
	Btex:bind(1)
	gradtex:bind(2)

	gl.glUniform1f(shader.uniforms.alpha.loc, 1)

	drawReading{
		geom = geom,
		lat = {42, 52, 45.021},
		lon = {74, 34, 16.004},
	}

	drawReading{
		geom = geom,
		lat = {33, 59, 38},		-- lat
		lon = {-80, -27, -56},	-- lon
	}

	if doDrawVectorField then
		drawVectorField(geom)
	end

	gl.glUniform1f(shader.uniforms.alpha.loc, drawAlpha[0])

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
	bool[0] = self.view.ortho
	if ig.igCheckbox('ortho', bool) then
		self.view.ortho = bool[0]
	end

	bool[0] = doDrawVectorField
	if ig.igCheckbox('draw vector field', bool) then
		doDrawVectorField = bool[0]
	end

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

	ig.igInputFloat('alpha', drawAlpha)
end

App():run()

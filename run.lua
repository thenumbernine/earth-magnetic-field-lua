#!/usr/bin/env luajit
local ffi = require 'ffi'
local vec3f = require 'vec-ffi.vec3f'
local vec4f = require 'vec-ffi.vec4f'
local quatf = require 'vec-ffi.quatf'
local template = require 'template'
local class = require 'ext.class'
local table = require 'ext.table'
local string = require 'ext.string'
local path = require 'ext.path'
local gl = require 'gl'
local glreport = require 'gl.report'
local ig = require 'imgui'
local GLTex2D = require 'gl.tex2d'
local GLProgram = require 'gl.program'
local glCallOrRun = require 'gl.call'
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
local year = 2020	-- TODO add support for dg/dh




-- load wmm
local lines = string.split(string.trim(path'wmm.cof':read()), '\n')
local header = string.split(string.trim(lines:remove(1)), '%s+')
assert(#header == 3)
print('model epoch', header[1])
print('model name', header[2])
print('date of release', header[3])
assert(lines:remove() == '999999999999999999999999999999999999999999999999')
assert(lines:remove() == '999999999999999999999999999999999999999999999999')

local wmm = {}
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



-- latitude = phi in [-pi/2, pi/2]
-- longitude = lambda in [-pi,pi]
-- returns xyz cartesian coordinates in units of the wgs84's 'a' parameter
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

-- expects xyz in cartesian units of wgs84's 'a' parameter
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
	-- TODO check domain and see if you can just use zp/r2D and sqrt(1-sinPhiSph^2)
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


local BStat = StatSet(
	'mag', 'x', 'y', 'z', 'mag2d',
	'div', 'div2d', 'curlZ', 'curlMag'
)

local londim -- dimension in 2D x dir / spherical phi / globe lambda dir
local latdim -- dimension in 2D y dir / spherical theta / globe phi dir

local Btex
local B2tex


local earthtex

local App = require 'imguiapp.withorbit'()

App.title = 'EM field'


local guivars = {
	geomIndex = 2,	-- north pole
	overlayIndex = 1,
	gradientIndex = 0,

	drawAlpha = 1,
	doDrawVectorField = true,

	fieldDT = 0,


	-- set to 0 to flatten vector field against surface
	fieldZScale = 1,

	arrowScale = 5,
}



local Geom = class()

function Geom:init(args)
	for k,v in pairs(args) do
		self[k] = v
	end
end

function Geom:draw()
	local height = 0
	local jres = 120
	local ires = 60
	--self.list = self.list or {}
	--glCallOrRun(self.list, function()
		for ibase=0,ires-1 do
			gl.glBegin(gl.GL_TRIANGLE_STRIP)
			for j=0,jres do
				local v = j/jres
				local lambda = math.rad((v * 2 - 1) * 180)
				for iofs=1,0,-1 do
					local i = ibase + iofs
					local u = i/ires
					local phi = math.rad((u * 2 - 1) * 90)
					
					local x,y,z = self:chart(phi, lambda, height)
					gl.glTexCoord2f(v, u)
					gl.glVertex3f(x, y, z)
				end
			end
			gl.glEnd()
		end
	--end)
end

local geoms = {
	Geom{
		name = 'Equirectangular',
		R = 2/math.pi,
		lambda0 = 0,
		phi0 = 0,
		phi1 = 0,
		updateGUI = function(self)
			ig.luatableInputFloat('R', self, 'R')
			ig.luatableInputFloat('lambda0', self, 'lambda0')
			ig.luatableInputFloat('phi0', self, 'phi0')
			ig.luatableInputFloat('phi1', self, 'phi1')
		end,
		chart = function(self, phi, lambda, height)
			return self.R * (lambda - self.lambda0) * math.cos(self.phi1),
				self.R * (phi - self.phi0),
				height / wgs84.a
		end,
		basis = function(self, phi, lambda, height)
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
	Geom{
		name = 'Azimuthal equidistant',
		chart = function(self, phi, lambda, height)
			local theta = .5 * math.pi - phi
			return
				math.cos(lambda) * theta,
				math.sin(lambda) * theta,
				height / wgs84.a
		end,
		basis = function(self, phi, lambda, height)
			local cosLambda = math.cos(lambda)
			local sinLambda = math.sin(lambda)
			return
				vec3f(-cosLambda, -sinLambda, 0),	-- d/dphi
				vec3f(-sinLambda, cosLambda, 0),	-- d/dlambda
				vec3f(0, 0, -1)						-- d/dheight
		end,
	},
	Geom{
		name = 'Mollweide',
		R = math.pi / 4,
		lambda0 = 0,	-- in degrees
		updateGUI = function(self)
			ig.luatableInputFloat('R', self, 'R')
			ig.luatableInputFloat('lambda0', self, 'lambda0')
		end,
		chart = function(self, phi, lambda, height)
			local theta = phi
			for i=1,10 do
				local dtheta = (2 * theta + math.sin(2 * theta) - math.pi * math.sin(phi)) / (2 + 2 * math.cos(theta))
				if math.abs(dtheta) < 1e-5 then break end
				theta = theta - dtheta
			end
			
			lambda = lambda - math.rad(self.lambda0)
			--[[ dumb and not working
			while lambda > math.pi do lambda = lambda - 2*math.pi end
			while lambda < -math.pi do lambda = lambda + 2*math.pi end
			--]]

			local x = self.R * math.sqrt(8) / math.pi * lambda * math.cos(theta)
			local y = self.R * math.sqrt(2) * math.sin(theta)
			return x, y, height / wgs84.a
		end,
		basis = function(self, phi, lambda, height)
			return
				vec3f(0, 1, 0),
				vec3f(1, 0, 0),
				vec3f(0, 0, -1)
		end,
	},
	Geom{
		name = 'WGS84',
		chart = function(self, phi, lambda, height)
			return latLonToCartesianWGS84(phi, lambda, height)
		end,
		basis = function(self, phi, lambda, height)
			return latLonToCartesianTangentSpaceWGS84(phi, lambda, height)
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



	local calc_b_shader = path'calc_b.shader':read()

	local vertexCode = [[
varying vec2 texcoordv;

void main() {
	texcoordv = gl_MultiTexCoord0.xy;
	gl_Position = ftransform();
}
]]
	local calcBFragmentCode = template(
[[
#version 120

uniform float dt;

]]..calc_b_shader..[[

#define M_PI <?=('%.49f'):format(math.pi)?>

varying vec2 texcoordv;

void main() {
	float phi = (texcoordv.y - .5) * M_PI;			//[-pi/2, pi/2]
	float lambda = (texcoordv.x - .5) * 2. * M_PI;	//[-pi, pi]
	vec3 B = calcB(vec3(phi, lambda, 0.));
	gl_FragColor = vec4(B, 1.);
}
]],
		{
			dt = 0,
			wgs84 = wgs84,
			wmm = wmm,
		}
	)

	path'calc_b.postproc.frag':write(calcBFragmentCode)





	do
		print'generating B field...'
		
		-- can be arbitrary
		-- but the WMM model is for 15 mins, so [360,180] x4
		londim = 1440
		latdim = 720
		
		local fbo = require 'gl.fbo'()
glreport'here'

		local calcBShader = GLProgram{
			vertexCode = vertexCode,
			fragmentCode = calcBFragmentCode,
		}
glreport'here'
		calcBShader:useNone()
glreport'here'

		Btex = GLTex2D{
			internalFormat = gl.GL_RGBA32F,
			width = londim,
			height = latdim,
			format = gl.GL_RGB,
			type = gl.GL_FLOAT,
			minFilter = gl.GL_LINEAR,
			magFilter = gl.GL_NEAREST,
			wrap = {
				s = gl.GL_REPEAT,
				t = gl.GL_REPEAT,
			},
		}
glreport'here'
	
		fbo:draw{
			viewport = {0, 0, londim, latdim},
			resetProjection = true,
			shader = calcBShader,
			dest = Btex,
		}
glreport'here'
	
		Btex:bind()
		Btex:generateMipmap()
		Btex:unbind()
glreport'here'

-- [=[ hmm, better way than copy paste?
		
		print'generating div B and curl B...'

		local calcB2Shader = GLProgram{
			vertexCode = vertexCode,
			fragmentCode = template(
[[
#version 120

uniform float dt;

]]..calc_b_shader..[[

#define M_PI <?=('%.49f'):format(math.pi)?>

#define latdim	<?=clnumber(latdim)?>
#define londim	<?=clnumber(londim)?>

vec3 dphi = vec3(M_PI / londim, 0., 0.);
vec3 dlambda = vec3(0., 2. * M_PI / latdim, 0.);
vec3 dheight = vec3(0., 0., 1.);

varying vec2 texcoordv;

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

	gl_FragColor = vec4(
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
					wgs84 = wgs84,
					wmm = wmm,
					dt = 0,
				}
			)
		}
glreport'here'
		calcB2Shader:useNone()
glreport'here'

		B2tex = GLTex2D{
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
		}
glreport'here'
	
		fbo:draw{
			viewport = {0, 0, londim, latdim},
			resetProjection = true,
			shader = calcB2Shader,
			dest = B2tex,
		}
glreport'here'


		-- only used for stat calc
		local Bdata = ffi.new('vec3f_t[?]', londim * latdim)
		Btex:toCPU(Bdata)
glreport'here'
		
glreport'here'
		
		--B2tex:generateMipmap()

		local B2data = ffi.new('vec4f_t[?]', londim * latdim)
		B2tex:bind()
		B2tex:toCPU(B2data)
glreport'here'
		
		B2tex:unbind()
glreport'here'

		local statgens = table{
			function(B, B2) return B:length() end,
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
				local B = Bdata[e]
				local B2 = B2data[e]
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
				local B = Bdata[e]
				local B2 = B2data[e]
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


	for _,overlay in ipairs(overlays) do
		overlay.shader = GLProgram{
			vertexCode = vertexCode,
			fragmentCode =
[[
#version 120
//120 needed for mat2x3 type

uniform float dt;

]]
..
template(
	calc_b_shader,
	{
		wgs84 = wgs84,
		wmm = wmm,
	}
)
..
template(
[[

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

varying vec2 texcoordv;

uniform sampler2D earthtex;
uniform sampler2D Btex;
uniform sampler2D B2tex;
uniform sampler1D hsvtex;
uniform float alpha;

mat3 latLonToCartesianTangentSpaceWGS84(vec3 plh) {
	float phi = plh.x;
	float lambda = plh.y;
	float height = plh.z;

	float cosLambda = cos(lambda);
	float sinLambda = sin(lambda);

	float cosPhi = cos(phi);
	float sinPhi = sin(phi);
	float dphi_cosPhi = -sinPhi;
	float dphi_sinPhi = cosPhi;

	float rCart = wgs84_a / sqrt(1. - wgs84_epssq * sinPhi * sinPhi);
	float tmp = sqrt(1. - wgs84_epssq * sinPhi * sinPhi);
	float dphi_rCart = wgs84_a / (tmp*tmp*tmp) * wgs84_epssq * sinPhi * dphi_sinPhi;

	float rCart_over_a = 1. / sqrt(1. - wgs84_epssq * sinPhi * sinPhi);

	float xp = (rCart + height) * cosPhi;
	float dphi_xp = dphi_rCart * cosPhi + (rCart + height) * dphi_cosPhi;
	float dheight_xp = cosPhi;

	float xp_over_a = (rCart_over_a + height / wgs84_a) * cosPhi;

	float zp = (rCart * (1. - wgs84_epssq) + height) * sinPhi;
	float dphi_zp = (dphi_rCart * (1. - wgs84_epssq)) * sinPhi + (rCart * (1. - wgs84_epssq) + height) * dphi_sinPhi;
	float dheight_zp = sinPhi;

	float zp_over_a = (rCart_over_a * (1. - wgs84_epssq) + height / wgs84_a) * sinPhi;

	float r2D = sqrt(xp * xp + zp * zp);
	float dphi_r2D = (xp * dphi_xp + zp * dphi_zp) / r2D;
	float dheight_r2D = (xp * dheight_xp + zp * dheight_zp) / r2D;

	float r2D_over_a = sqrt(xp_over_a * xp_over_a + zp_over_a * zp_over_a);
	float dphi_r2D_over_a = (xp_over_a * dphi_xp + zp_over_a * dphi_zp) / r2D;

	float sinPhiSph = zp / r2D;
	float dphi_sinPhiSph = (dphi_zp * r2D - zp * dphi_r2D) / (r2D * r2D);
	float dheight_sinPhiSph = (dheight_zp * r2D - zp * dheight_r2D) / (r2D * r2D);

	float cosPhiSph = sqrt(1. - sinPhiSph * sinPhiSph);
	// d/du sqrt(1 - x^2) = -x/sqrt(1 - x^2) dx/du;
	float dphi_cosPhiSph = -sinPhi / cosPhiSph * dphi_sinPhiSph;
	float dheight_cosPhiSph = -sinPhi / cosPhiSph * dheight_sinPhiSph;

	// float x = r2D * cosPhiSph / wgs84_a * cosLambda;
	// float y = r2D * cosPhiSph / wgs84_a * sinLambda;
	// float z = r2D * sinPhiSph / wgs84_a;

	float dphi_x = (dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * cosLambda;
	float dphi_y = (dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * sinLambda;
	float dphi_z = (dphi_r2D_over_a * sinPhiSph + r2D_over_a * dphi_sinPhiSph);

	float dlambda_x = -sinLambda;
	float dlambda_y = cosLambda;

	float dheight_x = (dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * cosLambda;
	float dheight_y = (dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * sinLambda;
	float dheight_z = (dheight_r2D * sinPhiSph + r2D * dheight_sinPhiSph);

	return mat3(
		vec3(dphi_x, dphi_y, dphi_z),
		vec3(dlambda_x, dlambda_y, 0),
		vec3(-dheight_x, -dheight_y, -dheight_z)
	);
}


void main() {
	float s = .5;
	float hsvBlend = .5;

	vec4 B2 = texture2D(B2tex, texcoordv);
#ifndef CALC_B_ON_GPU
	vec4 B = texture2D(Btex, texcoordv);
#else
	vec3 plh = vec3(
		(texcoordv.y - .5) * M_PI,			//phi
		(texcoordv.x - .5) * 2. * M_PI,		//lambda
		0.
	);
	vec3 B = calcB(plh);

#if 0
	mat3 e = latLonToCartesianTangentSpaceWGS84(plh);
	//TODO orthonormalize e
	//then TODO sample along each direction
	//then TODO convert from cartesian back to plh
#endif

#endif

	<?=overlay.code?>
	
	gl_FragColor = mix(
		texture2D(earthtex, vec2(texcoordv.x, 1. - texcoordv.y)),
		texture1D(hsvtex, s),
		hsvBlend);
	gl_FragColor.a = alpha;
}
]], 		{
				overlay = overlay,
				BStat = BStat,
				clnumber = clnumber,
			}),
			uniforms = {
				earthtex = 0,
				hsvtex = 1,
				Btex = 2,
				B2tex = 3,
				alpha = 1,
				dt = 0,
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

	gl.glEnable(gl.GL_CULL_FACE)
	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
end

local function degmintofrac(deg, min, sec)
	return deg + (1/60) * (min + (1/60) * sec)
end

local function drawReading(info)
	local geom = info.geom

	local height = 1e-2
	local lat = degmintofrac(table.unpack(info.lat))
	local lon = degmintofrac(table.unpack(info.lon))


	local phi = math.rad(lat)
	local lambda = math.rad(lon)
	local headingRad = math.rad(info.heading)

	local x,y,z = geom:chart(phi, lambda, height)
	local ex, ey, ez = geom:basis(phi, lambda, height)
	local H = (ex * math.cos(headingRad) + ey * math.sin(headingRad)) * .2
	
	gl.glColor3f(1,1,0)
	
	-- TODO geom should be transforms, and apply to all geom rendered
	gl.glPointSize(5)
	gl.glBegin(gl.GL_POINTS)
	
	gl.glVertex3f(x, y, z)
	
	gl.glEnd()
	gl.glPointSize(1)

	gl.glBegin(gl.GL_LINES)
	gl.glColor3f(0,0,1)
	gl.glVertex3f(x, y, z)
	gl.glVertex3f(x + ex.x * .2, y + ex.y * .2, z + ex.z * .2)
	gl.glColor3f(0,1,0)
	gl.glVertex3f(x, y, z)
	gl.glVertex3f(x + ey.x * .2, y + ey.y * .2, z + ey.z * .2)
	gl.glColor3f(1,0,1)
	gl.glVertex3f(x, y, z)
	gl.glVertex3f(x + H.x, y + H.y, z + H.z)
	gl.glEnd()
	
	-- now use wgs84 here regardless of 'geom'
	local pos = vec3f(latLonToCartesianWGS84(phi, lambda, height))
	local ex, ey, ez = latLonToCartesianTangentSpaceWGS84(phi, lambda, .1)
	local H = (ex * math.cos(headingRad) + ey * math.sin(headingRad)) * .2

	local axis = pos:cross(H):normalize()

	local degToSpan = 90
	local numSteps = 90
	local dtheta = degToSpan / numSteps
	local q = quatf():fromAngleAxis(axis.x, axis.y, axis.z, dtheta)
	
	-- draw north vs magnetic north heading
	gl.glColor3f(1,1,1)
	gl.glBegin(gl.GL_LINE_STRIP)

	gl.glVertex3f(geom:chart(cartesianToLatLonWGS84(pos:unpack())))
	for i=1,numSteps do
		
		q:rotate(pos, pos)
--print(pos, geom:chart(cartesianToLatLonWGS84(pos:unpack())))
		gl.glVertex3f(geom:chart(cartesianToLatLonWGS84(pos:unpack())))
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

local function drawVectorField(geom)

	local arrow = {
		{-.5, 0.},
		{.5, 0.},
		{.2, .3},
		{.5, 0.},
		{.2, -.3},
		{.5, 0.},
	}

	local height = 0
	local jres = 60
	local ires = 30
	local scale = guivars.arrowScale  / (BStat.mag.max * ires)
--	geom.list = geom.list or {}
--	glCallOrRun(geom.list, function()
	gl.glBegin(gl.GL_LINES)
	for i=0,ires do
		local u = i/ires
		local phi = math.rad((u * 2 - 1) * 90)
		for j=0,jres do
			local v = j/jres
			local lambda = math.rad((v * 2 - 1) * 180)
			local x,y,z = geom:chart(phi, lambda, 0)
			gl.glTexCoord2f(u,v)
			-- TODO calc B here, and add here
			local Bx, By, Bz = calcB(phi, lambda, height)
			Bx = Bx * scale	-- north component
			By = By * scale	-- east component
			Bz = Bz * scale	-- down component
			
			local ex, ey, ez = geom:basis(phi, lambda, height)

			--[[ verify the basis is correct
			scale = 3 / jres
			gl.glVertex3f(x, y, z)	gl.glVertex3f(x + ex.x * scale, y + ex.y * scale, z + ex.z * scale)
			gl.glVertex3f(x, y, z)	gl.glVertex3f(x + ey.x * scale, y + ey.y * scale, z + ey.z * scale)
			gl.glVertex3f(x, y, z)	gl.glVertex3f(x + ez.x * scale, y + ez.y * scale, z + ez.z * scale)
			--]]
			-- [[ draw our arrow
			local B = ex * Bx
				+ ey * By
				+ ez * Bz * guivars.fieldZScale

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


function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local shader = assert(overlays[tonumber(guivars.overlayIndex)+1]).shader
	local gradtex = assert(gradients[tonumber(guivars.gradientIndex)+1]).tex
	local geom = assert(geoms[tonumber(guivars.geomIndex)+1])

	shader:use()

	if shader.uniforms.dt then
		gl.glUniform1f(shader.uniforms.dt.loc, guivars.fieldDT)
	end

	earthtex:bind(0)
	gradtex:bind(1)
	Btex:bind(2)
	B2tex:bind(3)

	gl.glUniform1f(shader.uniforms.alpha.loc, 1)

	shader:useNone()

-- TODO more samples
--[[
	drawReading{
		geom = geom,
		lat = {42, 52, 45.021},
		lon = {74, 34, 16.004},
		heading = 5,
	}

	drawReading{
		geom = geom,
		lat = {33, 59, 38},		-- lat
		lon = {-80, -27, -56},	-- lon
		heading = -.5,
	}
--]]

	if guivars.doDrawVectorField then
		drawVectorField(geom)
	end
	
	shader:use()

	gl.glUniform1f(shader.uniforms.alpha.loc, guivars.drawAlpha)

	gl.glCullFace(gl.GL_FRONT)
	geom:draw()
	gl.glCullFace(gl.GL_BACK)
	geom:draw()

	B2tex:unbind(3)
	Btex:unbind(2)
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

	ig.igText'geom'
	for i,geom in ipairs(geoms) do
		ig.luatableRadioButton(geom.name, guivars, 'geomIndex', i-1)
	end

	ig.igSeparator()

	local geom = assert(geoms[tonumber(guivars.geomIndex)+1])
	if geom.updateGUI then
		geom:updateGUI()
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

	ig.luatableInputFloat('alpha', guivars, 'drawAlpha')
	ig.luatableInputFloat('field z', guivars, 'fieldZScale')
	ig.luatableInputFloat('field size', guivars, 'arrowScale')

	-- how linear are the g and h coeffs?
	-- can I just factor out the dt?
	--ig.luatableInputFloat('time from '..wmm.epoch, guivars, 'fieldDT')
	ig.luatableSliderFloat('time from '..tostring(wmm.epoch), guivars, 'fieldDT', 0, 5)
end

return App():run()

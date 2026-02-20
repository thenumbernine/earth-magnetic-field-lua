-- wmm loading stuff used by run.lua and build-btex.lua
local class = require 'ext.class'
local assert = require 'ext.assert'
local path = require 'ext.path'
local string = require 'ext.string'

local WMM = class()

--[[
args:
	cof = -cof.txt file, default 'wmm2025-cof.txt'
	nMax = nMax, default #wmm, which is 12 for wmm2025.cof
--]]
function WMM:init(args)

	-- Sets WGS-84 parameters
	local wgs84 = {}
	wgs84.a = 6378.137 --semi-major axis of the ellipsoid in
	wgs84.b = 6356.7523142 --semi-minor axis of the ellipsoid in
	wgs84.fla = 1 / 298.257223563 -- flattening
	wgs84.eps = math.sqrt(1 - (wgs84.b * wgs84.b) / (wgs84.a * wgs84.a)) --first eccentricity
	wgs84.epssq = wgs84.eps * wgs84.eps --first eccentricity squared
	wgs84.re = 6371.2 -- Earth's radius
	self.wgs84 = wgs84


	--local HeightAboveEllipsoid = 0		-- compute at z=0 for now
	-- TODO add support for dg/dh

	local wmmfn = args.cof or 'wmm2025-cof.txt'	-- 'wmm2020-cof.txt'

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

	self.nMax = args.nMax or #wmm
	assert.le(self.nMax, #wmm, "nMax must be <= #wmm")
	self.wmm = wmm
end

-- ported from WMM2020 GeomagnetismLibrary.c
-- phi = radians
-- lambda = radians
-- height = meters
function WMM:calcB(phi, lambda, height)
	local wmm = self.wmm
	local wgs84 = self.wgs84
	local nMax = self.nMax

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
function WMM:cartesianToLatLonWGS84(x, y, z)
	local wgs84 = self.wgs84

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

return WMM

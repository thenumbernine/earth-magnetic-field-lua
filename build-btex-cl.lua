#!/usr/bin/env luajit
--[[
The calc_b.shader takes my 10yo nvidia computer about 1 minute to link when using nMax == #wmm == 12
It will probably load a high resolution FITS file of the B field much quicker.
So let's just save one here.

I did it lazy.
- you could speed it up by using cached neighbors
- you could speed it up by multithreading it
- you could speed it up by writing it as a compute-gpu kernel (tho the nvidia linker is taking 1minute to link a glsl program made for nMax=#wmm=12. ....)

but then I did all those things to try to fix the B gradient error at the poles
TODO FIX THE B GRADIENT ERROR AT THE POLES

TODO this in opencl because lua is slow
--]]
local ffi = require 'ffi'
local assert = require 'ext.assert'
local template = require 'template'
local math = require 'ext.math'
local path = require 'ext.path'
local cmdline = require 'ext.cmdline'(...)
local WMM = require 'earth-magnetic-field.wmm'
local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'

local W = WMM{
	cof = cmdline.cof,
	nMax = cmdline.nMax,
}
local wgs84 = W.wgs84

local Image = require 'image'

local londim = 4096
local latdim = 2048

local env = CLEnv{
	useGLSharing = false,
	getPlatform = CLEnv.getPlatformFromCmdLine(...),
	getDevices = CLEnv.getDevicesFromCmdLine(...),
	deviceType = CLEnv.getDeviceTypeFromCmdLine(...),
	size = {londim, latdim},
	real = 'double',
}

local dphi = math.pi / londim
local dlambda = 2 * math.pi / latdim
local dheight = 1000

local code = template([[

<?
nMax = nMax or #wmm
local glnumber = require 'gl.number'
local function int(x) return x < 0 and math.ceil(x) or math.floor(x) end
for k,v in pairs(wgs84) do
?>#define wgs84_<?=k?>	<?=glnumber(v)?>
<?
end
?>
#define nMax		<?=int(nMax)?>
#define numTerms 	((nMax + 1) * (nMax + 2) / 2)

// why is OpenCL like this?
#define _real2(a,b) (real2)(a,b)
#define _real4(a,b,c,d) (real4)(a,b,c,d)

real2 cplxmul(real2 a, real2 b) {
	return _real2(
		a.x * b.x - a.y * b.y,
		a.x * b.y + a.y * b.x
	);
}

/*
this is called "shit-multiply" because the folks at OpenGL named their matrixes backwards.
a is the diagonal of a 3x3 matrix,
b is a 3x2 matrix, but GLSL calls it a 2x2.
c is a 2x1 vector
*/
real4 glslStyle_vec3_mat2x3_vec2_shitMul(
	real4 a,
	real4 bcol1,
	real4 bcol2,
	real2 c
) {
	real4 r = bcol1 * c.x + bcol2 * c.y;
	// opencl per elem mul?
	return _real4(a.x * r.x, a.y * r.y, a.z * r.z, 0.);
}

// ported from WMM2020 GeomagnetismLibrary.c
// plh = phi, lambda, height in meters, dt in some time unit
// phi in [-π/2, π/2]
// lambda in [-π, π]
// height is in m
real4 calcB(real4 plh) {
	plh.z *= 1e-3;	// m to km
	real dt = plh.w;

	// begin MAG_GeodeticToSpherical
	real2 cisPhi = _real2(cos(plh.x), sin(plh.x));

	// convert from geodetic WGS-84 to spherical coordiantes
	real rc = wgs84_a / sqrt(1. - wgs84_epssq * cisPhi.y * cisPhi.y);

	real2 xzp = _real2(
		rc + plh.z,
		rc * (1. - wgs84_epssq) + plh.z
	) * cisPhi;

	// spherical results:

	real invR = 1. / length(xzp);

	real2 cisPhiSph;
	cisPhiSph.y = xzp.y * invR;	// geocentric latitude sin & cos
	cisPhiSph.x = sqrt(1. - cisPhiSph.y * cisPhiSph.y);

	// longitude is the same
	// end MAG_GeodeticToSpherical

	// begin MAG_Geomag
	// begin MAG_ComputeSphericalHarmonicVariables

	real2 cisLambda = _real2(cos(plh.y), sin(plh.y));

	//real2 cisLambdaToTheM[nMax+1];
	real2 cisLambdaToTheM_0 = _real2(1., 0.);
	real2 cisLambdaToTheM_1 = cisLambda;
<? for m=2,nMax do
?>	real2 cisLambdaToTheM_<?=m?> = cplxmul(cisLambdaToTheM_<?=m-1?>, cisLambda);
<? end

	-- end MAG_ComputeSphericalHarmonicVariables
	-- begin MAG_AssociatedLegendreFunction

	-- Compute the ration between the the Schmidt quasi-normalized associated Legendre
	-- functions and the Gauss-normalized version.

	local schmidtQuasiNorm = {}
	schmidtQuasiNorm[0] = 1

	for n=1,nMax do
		local index = (n * (n + 1) / 2)
		local index1 = (n - 1) * n / 2

		-- for m = 0
		schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (2 * n - 1) / n

		for m=1,n do
			local index = (n * (n + 1) / 2 + m)
			local index1 = (n * (n + 1) / 2 + m - 1)
			schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * math.sqrt( ((n - m + 1) * (m == 1 and 2 or 1)) / (n + m))
		end
	end

	-- begin MAG_PcupLow
?>

	//real2 P[numTerms];	//Legendre function & derivative
	real2 P_0 = _real2(1., 0.);
<?
	--	 First,	Compute the Gauss-normalized associated Legendre functions

	for n=1,nMax do
		for m=0,n do
			local index = n * (n + 1) / 2 + m
			if n == m then
				local index1 = (n - 1) * n / 2 + m - 1
?>	real2 P_<?=int(index)?> = _real2(cisPhiSph.x * P_<?=int(index1)?>.x, cisPhiSph.x * P_<?=int(index1)?>.y + cisPhiSph.y * P_<?=int(index1)?>.x);
<?
			elseif n == 1 and m == 0 then
				local index1 = (n - 1) * n / 2 + m
?>	real2 P_<?=int(index)?> = _real2(cisPhiSph.y * P_<?=int(index1)?>.x, cisPhiSph.y * P_<?=int(index1)?>.y - cisPhiSph.x * P_<?=int(index1)?>.x);
<?
			elseif n > 1 and n ~= m then
				local index1 = (n - 2) * (n - 1) / 2 + m
				local index2 = (n - 1) * n / 2 + m
				if m > n - 2 then
?>	real2 P_<?=int(index)?> = _real2(cisPhiSph.y * P_<?=int(index2)?>.x, cisPhiSph.y * P_<?=int(index2)?>.y - cisPhiSph.x * P_<?=int(index2)?>.x);
<?
				else
					local k = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3))
?>	real2 P_<?=int(index)?> = _real2(cisPhiSph.y * P_<?=int(index2)?>.x - <?=glnumber(k)?> * P_<?=int(index1)?>.x, cisPhiSph.y * P_<?=int(index2)?>.y - cisPhiSph.x * P_<?=int(index2)?>.x - <?=glnumber(k)?> * P_<?=int(index1)?>.y);
<?
				end
			end
		end
	end

-- this is now baked into the B sum:

	-- Converts the Gauss-normalized associated Legendre
	-- functions to the Schmidt quasi-normalized version using pre-computed
	-- relation stored in the variable schmidtQuasiNorm

	-- The sign is changed since the new WMM routines use derivative with respect to latitude
	-- insted of co-latitude

	-- end MAG_PcupLow
	-- end MAG_AssociatedLegendreFunction
	-- begin MAG_Summation

?>

	real earthRadOverR = wgs84_re * invR;

	real4 B = earthRadOverR * earthRadOverR * (
		_real4(0., 0., 0., 0.)
<?
	for n=1,nMax do
?>		+ earthRadOverR * (
			_real4(0., 0., 0., 0.)
<?
		for m=0,n do
			local index = (n * (n + 1) / 2 + m)

			--.g .h .gt .ht == .xyzw
			-- then again, looks like I'm not doing any gt/ht calculations...
			-- that means my reading is strictly 2020, right?

			--		   nMax  (n+2) n     m            m           m
			--		Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1         m=0   n            n           n
			-- Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius.

			--		  1 nMax  (n+2)    n     m            m           m
			--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1             m=0   n            n           n
			-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.

			--		    nMax  	(n+2) 	  n     m            m           m
			--		Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
			--						n=1      	      m=0   n            n           n
			-- Equation 12 in the WMM Technical report.  Derivative with respect to radius.

?> 			+ glslStyle_vec3_mat2x3_vec2_shitMul(
<?			-- rhs col = lhs row mul ...
			if m == 0 then
?> 				_real4(-P_<?=int(index)?>.y, 0., -P_<?=int(index)?>.x, 0.),
<?			else
?> 				_real4(-P_<?=int(index)?>.y, P_<?=int(index)?>.x, -P_<?=int(index)?>.x, 0.),
<?			end 
		
			-- 3x2 matrix col #1:
?>				_real4(
					<?=glnumber(wmm[n][m].g * -schmidtQuasiNorm[index])?>
					- dt * <?=glnumber(wmm[n][m].dg_dt * -schmidtQuasiNorm[index])?>,
<? 			if m == 0 then 
?>					0.,
<? 			else 
?>					<?=glnumber(-wmm[n][m].h * m * schmidtQuasiNorm[index])?>
					- dt * <?=glnumber(wmm[n][m].dh_dt * m * schmidtQuasiNorm[index])?>,
<? 			end 
?> 					<?=glnumber(wmm[n][m].g * (n + 1) * schmidtQuasiNorm[index])?>
					- dt * <?=glnumber(wmm[n][m].dg_dt * (n + 1) * schmidtQuasiNorm[index])?>,
					0.
				), 
<? 			-- 3x2 matrix col #2:
?>				_real4(
					<?=glnumber(wmm[n][m].h * -schmidtQuasiNorm[index])?>
					- dt * <?=glnumber(wmm[n][m].dh_dt * -schmidtQuasiNorm[index])?>,
<? 			if m == 0 then 
?>					0.,
<?			else 
?>					<?=glnumber(wmm[n][m].g * m * schmidtQuasiNorm[index])?>
					- dt * <?=glnumber(wmm[n][m].dg_dt * m * schmidtQuasiNorm[index])?>,
<? 			end 
?>
					<?=glnumber(wmm[n][m].h * (n + 1) * schmidtQuasiNorm[index])?>
					- dt * <?=glnumber(wmm[n][m].dh_dt * (n + 1) * schmidtQuasiNorm[index])?>,
					0.
				),
<?			-- right mul	
?>			 	cisLambdaToTheM_<?=int(m)?>
			)
<?		end
	end
?>
	<?for n=0,nMax do
?>)<?
	end
?>;

	if (cisPhiSph.x < -1e-10 || cisPhiSph.x > 1e-10) {
		B.y /= cisPhiSph.x;
	} else {
		// Special calculation for component - By - at Geographic poles.
		// If the user wants to avoid using this function, please make sure that
		// the latitude is not exactly +/-90. An option is to make use the function
		// MAG_CheckGeographicPoles.
		// begin MAG_SummationSpecial

		//real PS[nMax+1];
		real PS_0 = 1.;

		B.y = 0.;

		real earthRadOverRToTheN = earthRadOverR * earthRadOverR;
<?
		local schmidtQuasiNorm1 = 1
		for n=1,nMax do
?>		earthRadOverRToTheN *= earthRadOverR;
<?
			--Compute the ration between the Gauss-normalized associated Legendre
			-- functions and the Schmidt quasi-normalized version. This is equivalent to
			-- sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)!

			local m = 1
			local schmidtQuasiNorm2 = schmidtQuasiNorm1 * (2 * n - 1) / n
			local schmidtQuasiNorm3 = schmidtQuasiNorm2 * math.sqrt((n * 2) / (n + 1))
			local schmidtQuasiNorm1 = schmidtQuasiNorm2
			if n == 1 then
?>		real PS_<?=int(n)?> = PS_<?=int(n-1)?>;
<? 			else
				local k = (((n - 1) * (n - 1)) - 1) / ((2 * n - 1) * (2 * n - 3))
?>		real PS_<?=int(n)?> = cisPhiSph.y * PS_<?=int(n-1)?> - <?=glnumber(k)?> * PS_<?=int(n-2)?>;
<?
			end

			--		  1 nMax  (n+2)    n     m            m           m
			--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1             m=0   n            n           n
			-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.

?>		B.y += earthRadOverRToTheN * dot(cisLambda, _real2(<?=
			glnumber(-wmm[n][m].h * schmidtQuasiNorm3)
		?>, <?=
			glnumber(wmm[n][m].g * schmidtQuasiNorm3)
		?>)) * PS_<?=int(n)?>;
<?
		end

		-- end MAG_SummationSpecial
?>	}

	// end MAG_Summation
	// end MAG_Geomag

	return B;

}

kernel void calcBBuf(
	global float4 * BBuf
) {
	initKernel();

	const size_t latdim = size.y;
	real v = ((real)i.y + .5) / (real)latdim;
	real lat = (v * 2. - 1.) * 90.;
	real phi = lat * M_PI / 180.;
	real cos_phi = cos(phi);

	const size_t londim = size.x;
	real u = ((real)i.x + .5) / (real)londim;
	real lon = (u * 2. - 1.) * 180.;
	real lambda = lon * M_PI / 180.;

	BBuf[index] = convert_float4(calcB(_real4(phi, lambda, 0., 0.)));
}
]],	{
		wgs84 = W.wgs84,
		wmm = W.wmm,
		nMax = W.nMax,
	})
path'build-btex.cl':write(code)
local program = env:program{code = code}
program:compile()

local BBuf = env:buffer{name='BBuf', type='float4'}
local calcBKernel = program:kernel'calcBBuf'
calcBKernel.obj:setArg(0, BBuf.obj)
calcBKernel()
local BImg = Image(londim, latdim, 4, 'float')	-- Bx, By, Bz, 0
BBuf:toCPU(BImg.buffer)
BImg:save'B.fits'
error'FINSIHME'

local B2Buf = env:buffer{name='B2Buf', type='real4'}
local B3Buf = env:buffer{name='B3Buf', type='real4'}

print'Building B'
local lastTime = os.time()
local Bptr = BImg.buffer+0
for j=0,latdim-1 do
	for i=0,londim-1 do
		local u = (i + .5) / londim
		local lon = (u * 2 - 1) * 180
		local lambda = math.rad(lon)

		local thisTime = os.time()
		if lastTime ~= thisTime then
			lastTime = thisTime
			print(
				(100 * (i + londim * j) / (latdim * londim))
				..'%% complete'
			)
		end

		-- calcB

		local Bx, By, Bz = W:calcB(phi, lambda, 0)

		Bptr[0] = Bx	Bptr=Bptr+1
		Bptr[0] = By	Bptr=Bptr+1
		Bptr[0] = Bz	Bptr=Bptr+1
		Bptr[0] = 0		Bptr=Bptr+1
	end
end
assert.eq(Bptr, BImg.buffer + BImg.channels * BImg.width * BImg.height)


print'Building B2 and B3'
local B2Img = Image(londim, latdim, 4, 'float')	-- div B, div2D B, curl B, |curl B|
local B3Img = Image(londim, latdim, 4, 'float')	-- |B| |B|2d, 0, 0

local lastTime = os.time()
local Bptr = BImg.buffer+0
local B2ptr = B2Img.buffer+0
local B3ptr = B3Img.buffer+0
for j=0,latdim-1 do
	local v = (j + .5) / latdim
	local lat = (v * 2 - 1) * 90
	local phi = math.rad(lat)
	local cos_phi = math.cos(phi)
	for i=0,londim-1 do
		local u = (i + .5) / londim
		local lon = (u * 2 - 1) * 180
		local lambda = math.rad(lon)

		local thisTime = os.time()
		if lastTime ~= thisTime then
			lastTime = thisTime
			print(
				(100 * (i + londim * j) / (latdim * londim))
				..'%% complete'
			)
		end

		local Bx, By, Bz = Bptr[0], Bptr[1], Bptr[2]

		-- calcB2

		-- [[ calculated again
		local Bx_phiR, By_phiR, Bz_phiR = W:calcB(phi + dphi, lambda, 0)
		local Bx_phiL, By_phiL, Bz_phiL = W:calcB(phi - dphi, lambda, 0)

		local Bx_lambdaR, By_lambdaR, Bz_lambdaR = W:calcB(phi, lambda + dlambda, 0)
		local Bx_lambdaL, By_lambdaL, Bz_lambdaL = W:calcB(phi, lambda - dlambda, 0)
		--]]
		--[[ use cached BImg copy
		local flip
		local iR = i
		local jR = j + 1
		if jR > latdim-1 then
			-- across the pole?
			-- lat -> lat + 180 degrees
			-- lon -> -lon
			-- flip everything
			iR = (iR + bit.rshift(londim, 1)) % londim
			jR = latdim-1 - (jR - (latdim-1))
			flip = true
		end
assert.le(0, iR) assert.lt(iR, londim)
assert.le(0, jR) assert.lt(jR, latdim)
		local Bx_phiR = BImg.buffer[0 + BImg.channels * (iR + londim * jR)]
		local By_phiR = BImg.buffer[1 + BImg.channels * (iR + londim * jR)]
		local Bz_phiR = BImg.buffer[2 + BImg.channels * (iR + londim * jR)]
		if flip then -- I flip these when I go over the pole ... right?
			Bx_phiR = -Bx_phiR
			By_phiR = -By_phiR
		end

		local iL = i
		local jL = j - 1
		local flip
		if jL < 0 then
			-- same
			iL = (iL + bit.rshift(londim, 1)) % londim
			jL = -jL
			flip = true
		end
assert.le(0, iL) assert.lt(iL, londim)
assert.le(0, jL) assert.lt(jL, latdim)
		local Bx_phiL = BImg.buffer[0 + BImg.channels * (iL + londim * jL)]
		local By_phiL = BImg.buffer[1 + BImg.channels * (iL + londim * jL)]
		local Bz_phiL = BImg.buffer[2 + BImg.channels * (iL + londim * jL)]
		if flip then -- I flip these when I go over the pole ... right?
			Bx_phiL = -Bx_phiL
			By_phiL = -By_phiL
		end

		local jR = j
		local iR = i + 1
		if iR > londim-1 then
			iR = iR % londim
		end
assert.le(0, iR) assert.lt(iR, londim)
assert.le(0, jR) assert.lt(jR, latdim)
		local Bx_lambdaR = BImg.buffer[0 + BImg.channels * (iR + londim * jR)]
		local By_lambdaR = BImg.buffer[1 + BImg.channels * (iR + londim * jR)]
		local Bz_lambdaR = BImg.buffer[2 + BImg.channels * (iR + londim * jR)]

		local jL = j
		local iL = i - 1
		if iL < 0 then
			iL = iL % londim
		end
assert.le(0, iL) assert.lt(iL, londim)
assert.le(0, jL) assert.lt(jL, latdim)
		local Bx_lambdaL = BImg.buffer[0 + BImg.channels * (iL + londim * jL)]
		local By_lambdaL = BImg.buffer[1 + BImg.channels * (iL + londim * jL)]
		local Bz_lambdaL = BImg.buffer[2 + BImg.channels * (iL + londim * jL)]
		--]]

		-- there is no separate altitude cache
		local Bx_heightR, By_heightR, Bz_heightR = W:calcB(phi, lambda, 0 + dheight)
		local Bx_heightL, By_heightL, Bz_heightL = W:calcB(phi, lambda, 0 - dheight)

		-- [[ TODO both here and run.lua's calcB2Code, I am probably off in my magnitudes because it looks to be all zeroes
		local dBx_dphi = (Bx_phiR - Bx_phiL) / (2 * dphi) / (wgs84.a * 1e+3 * cos_phi)
		local dBy_dphi = (By_phiR - By_phiL) / (2 * dphi) / (wgs84.a * 1e+3 * cos_phi)
		local dBz_dphi = (Bz_phiR - Bz_phiL) / (2 * dphi) / (wgs84.a * 1e+3 * cos_phi)

		local dBx_dlambda = (Bx_lambdaR - Bx_lambdaL) / (2 * dlambda) / (wgs84.a * 1e+3)
		local dBy_dlambda = (By_lambdaR - By_lambdaL) / (2 * dlambda) / (wgs84.a * 1e+3)
		local dBz_dlambda = (Bz_lambdaR - Bz_lambdaL) / (2 * dlambda) / (wgs84.a * 1e+3)

		local dBx_dheight = (Bx_heightR - Bx_heightL) / (2 * dheight)
		local dBy_dheight = (By_heightR - By_heightL) / (2 * dheight)
		local dBz_dheight = (Bz_heightR - Bz_heightL) / (2 * dheight)
		--]]

		local div2D_B = dBx_dphi + dBy_dlambda
		local div_B = div2D_B + dBz_dheight

		local curl_B_x = dBz_dlambda - dBy_dheight
		local curl_B_y = dBx_dheight - dBz_dphi
		local curl_B_z = dBy_dphi - dBx_dlambda

		local len_curl_B = math.sqrt(curl_B_x^2 + curl_B_y^2 + curl_B_z^2)

		B2ptr[0] = div_B		B2ptr=B2ptr+1
		B2ptr[0] = div2D_B		B2ptr=B2ptr+1
		B2ptr[0] = curl_B_z		B2ptr=B2ptr+1
		B2ptr[0] = len_curl_B	B2ptr=B2ptr+1


		-- calcB3


		local len_B2 = math.sqrt(Bx^2 + By^2)
		local len_B3 = math.sqrt(len_B2^2 + Bz^2)

		B3ptr[0] = len_B2		B3ptr=B3ptr+1
		B3ptr[0] = len_B3		B3ptr=B3ptr+1
		B3ptr[0] = 0			B3ptr=B3ptr+1
		B3ptr[0] = 0			B3ptr=B3ptr+1


		Bptr=Bptr+4
	end
end
assert.eq(Bptr, BImg.buffer + BImg.channels * BImg.width * BImg.height)
assert.eq(B2ptr, B2Img.buffer + B2Img.channels * B2Img.width * B2Img.height)
assert.eq(B3ptr, B3Img.buffer + B3Img.channels * B3Img.width * B3Img.height)

B2Img:save'B2.fits'
B3Img:save'B3.fits'

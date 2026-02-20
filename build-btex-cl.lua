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

kernel void calcBBufs(
	global float4 * BBuf,
	global float4 * B2Buf,
	global float4 * B3Buf
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

	real4 plh = _real4(phi, lambda, 0., 0.);

	real4 B = calcB(plh);
	BBuf[index] = convert_float4(B);


	real4 dphi = _real4(M_PI / (real)londim, 0., 0., 0.);
	real4 dlambda = _real4(0., 2. * M_PI / (real)latdim, 0., 0.);
	real4 dheight = _real4(0., 0., 1000., 0.);

	// TODO units anyone?
	real4 dphi_B = (calcB(plh + dphi) - calcB(plh - dphi)) / (2. * dphi.x) / (wgs84_a * 1e+3 * cos(plh.x));
	real4 dlambda_B = (calcB(plh + dlambda) - calcB(plh - dlambda)) / (2. * dlambda.y) / (wgs84_a * 1e+3);
	real4 dheight_B = (calcB(plh + dheight) - calcB(plh - dheight)) / (2. * dheight.z);

	real div2D_B = dphi_B.x + dlambda_B.y;
	real div_B = div2D_B + dheight_B.z;

	real4 curl_B = _real4(
		dlambda_B.z - dheight_B.y,
		dheight_B.x - dphi_B.z,
		dphi_B.y - dlambda_B.x,
		0.
	);

	B2Buf[index] = convert_float4(_real4(
		div_B,
		div2D_B,
		curl_B.z,
		length(curl_B)
	));
	
	B3Buf[index] = convert_float4(_real4(
		length(B.xyz),
		length(B.xy),
		0.,
		1.
	));
}
]],	{
		wgs84 = W.wgs84,
		wmm = W.wmm,
		nMax = W.nMax,
	})

path'build-btex.cl':write(code)
local program = env:program{code = code}
program:compile()
local calcBKernel = program:kernel'calcBBufs'

local BBuf = env:buffer{name='BBuf', type='float4'}
local B2Buf = env:buffer{name='B2Buf', type='float4'}
local B3Buf = env:buffer{name='B3Buf', type='float4'}

calcBKernel.obj:setArg(0, BBuf.obj)
calcBKernel.obj:setArg(1, B2Buf.obj)
calcBKernel.obj:setArg(2, B3Buf.obj)
calcBKernel()

local BImg = Image(londim, latdim, 4, 'float')	-- Bx, By, Bz, 0
local B2Img = Image(londim, latdim, 4, 'float')	-- div B, div2D B, curl B, |curl B|
local B3Img = Image(londim, latdim, 4, 'float')	-- |B| |B|2d, 0, 0

BBuf:toCPU(BImg.buffer)
B2Buf:toCPU(B2Img.buffer)
B3Buf:toCPU(B3Img.buffer)

BImg:save'B.fits'
B2Img:save'B2.fits'
B3Img:save'B3.fits'

BImg:normalize():rgb():save'B.png'
B2Img:normalize():rgb():save'B2.png'
B3Img:normalize():rgb():save'B3.png'

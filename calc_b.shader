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

vec2 cplxmul(vec2 a, vec2 b) {
	return vec2(
		a.x * b.x - a.y * b.y,
		a.x * b.y + a.y * b.x
	);
}

/*
plh.xyz = phi, lambda, height

phi in [-pi/2, pi/2]
lambda in [-pi, pi]
height is in m
*/
vec3 calcB(vec3 plh) {
	plh.z *= 1e-3;	// m to km

	// begin MAG_GeodeticToSpherical
	vec2 cisPhi = vec2(cos(plh.x), sin(plh.x));

	// convert from geodetic WGS-84 to spherical coordiantes
	float rc = wgs84_a / sqrt(1. - wgs84_epssq * cisPhi.y * cisPhi.y);

	vec2 xzp = vec2(
		rc + plh.z,
		rc * (1. - wgs84_epssq) + plh.z
	) * cisPhi;

	// spherical results:

	float invR = 1. / length(xzp);

	vec2 cisPhiSph;
	cisPhiSph.y = xzp.y * invR;	// geocentric latitude sin & cos
	cisPhiSph.x = sqrt(1. - cisPhiSph.y * cisPhiSph.y);

	// longitude is the same
	// end MAG_GeodeticToSpherical

	// begin MAG_Geomag
	// begin MAG_ComputeSphericalHarmonicVariables

	vec2 cisLambda = vec2(cos(plh.y), sin(plh.y));

	//vec2 cisLambdaToTheM[nMax+1];
	vec2 cisLambdaToTheM_0 = vec2(1., 0.);
	vec2 cisLambdaToTheM_1 = cisLambda;
<? for m=2,nMax do
?>	vec2 cisLambdaToTheM_<?=m?> = cplxmul(cisLambdaToTheM_<?=m-1?>, cisLambda);
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

	//vec2 P[numTerms];	//Legendre function & derivative
	vec2 P_0 = vec2(1., 0.);
<?
	--	 First,	Compute the Gauss-normalized associated Legendre functions

	for n=1,nMax do
		for m=0,n do
			local index = n * (n + 1) / 2 + m
			if n == m then
				local index1 = (n - 1) * n / 2 + m - 1
?>	vec2 P_<?=int(index)?> = vec2(cisPhiSph.x * P_<?=int(index1)?>.x, cisPhiSph.x * P_<?=int(index1)?>.y + cisPhiSph.y * P_<?=int(index1)?>.x);
<?
			elseif n == 1 and m == 0 then
				local index1 = (n - 1) * n / 2 + m
?>	vec2 P_<?=int(index)?> = vec2(cisPhiSph.y * P_<?=int(index1)?>.x, cisPhiSph.y * P_<?=int(index1)?>.y - cisPhiSph.x * P_<?=int(index1)?>.x);
<?
			elseif n > 1 and n ~= m then
				local index1 = (n - 2) * (n - 1) / 2 + m
				local index2 = (n - 1) * n / 2 + m
				if m > n - 2 then
?>	vec2 P_<?=int(index)?> = vec2(cisPhiSph.y * P_<?=int(index2)?>.x, cisPhiSph.y * P_<?=int(index2)?>.y - cisPhiSph.x * P_<?=int(index2)?>.x);
<?
				else
					local k = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3))
?>	vec2 P_<?=int(index)?> = vec2(cisPhiSph.y * P_<?=int(index2)?>.x - <?=glnumber(k)?> * P_<?=int(index1)?>.x, cisPhiSph.y * P_<?=int(index2)?>.y - cisPhiSph.x * P_<?=int(index2)?>.x - <?=glnumber(k)?> * P_<?=int(index1)?>.y);
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

	float earthRadOverR = wgs84_re * invR;

	vec3 B = earthRadOverR * earthRadOverR * (vec3(0., 0., 0.)
<?
	for n=1,nMax do
?>		+ earthRadOverR * (
			vec3(0., 0., 0.)
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

?>			+ vec3(-P_<?=int(index)?>.y, <?
			if m == 0 then
				?>0.<?
			else
				?>P_<?=int(index)?>.x<? end ?>, -P_<?=int(index)?>.x) * (mat2x3(vec3(<?=
						glnumber(wmm[n][m].g * -schmidtQuasiNorm[index])
					?> - dt * <?=
						glnumber(wmm[n][m].dg_dt * -schmidtQuasiNorm[index])
					?>, <? if m == 0 then ?>0.<? else ?><?=
						glnumber(-wmm[n][m].h * m * schmidtQuasiNorm[index])
					?> - dt * <?=
						glnumber(wmm[n][m].dh_dt * m * schmidtQuasiNorm[index])
					?><? end ?>, <?=
						glnumber(wmm[n][m].g * (n + 1) * schmidtQuasiNorm[index])
					?> - dt * <?=
						glnumber(wmm[n][m].dg_dt * (n + 1) * schmidtQuasiNorm[index])
					?>), vec3(<?=
						glnumber(wmm[n][m].h * -schmidtQuasiNorm[index])
					?> - dt * <?=
						glnumber(wmm[n][m].dh_dt * -schmidtQuasiNorm[index])
					?>, <? if m == 0 then ?>0.<? else ?><?=
						glnumber(wmm[n][m].g * m * schmidtQuasiNorm[index])
					?> - dt * <?=
						glnumber(wmm[n][m].dg_dt * m * schmidtQuasiNorm[index])
					?><? end ?>, <?=
						glnumber(wmm[n][m].h * (n + 1) * schmidtQuasiNorm[index])
					?> - dt * <?=
						glnumber(wmm[n][m].dh_dt * (n + 1) * schmidtQuasiNorm[index])
					?>)) * cisLambdaToTheM_<?=int(m)?>)
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

		//float PS[numTerms];
		float PS_0 = 1.;

		B.y = 0.;

		float earthRadOverRToTheN = earthRadOverR * earthRadOverR;
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
?>		float PS_<?=int(n)?> = PS_<?=int(n-1)?>;
<? 			else
				local k = (((n - 1) * (n - 1)) - 1) / ((2 * n - 1) * (2 * n - 3))
?>		float PS_<?=int(n)?> = cisPhiSph.y * PS_<?=int(n-1)?> - <?=glnumber(k)?> * PS_<?=int(n-2)?>;
<?
			end

			--		  1 nMax  (n+2)    n     m            m           m
			--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1             m=0   n            n           n
			-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.

?>		B.y += earthRadOverRToTheN * dot(cisLambda, vec2(<?=
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

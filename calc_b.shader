<? 
local nMax = #wmm
local clnumber = require 'cl.obj.number' 
for k,v in pairs(wgs84) do
?>#define wgs84_<?=k?>	<?=clnumber(v)?>
<?
end
?>
#define nMax		<?=nMax?>
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
*/
vec3 calcB(vec3 plh) {
		
	// begin MAG_GeodeticToSpherical
	vec2 cisPhi = vec2(cos(plh.x), sin(plh.x));

	// convert from geodetic WGS-84 to spherical coordiantes
	float rc = wgs84_a / sqrt(1. - wgs84_epssq * cisPhi.y * cisPhi.y);
	
	vec2 xzp = vec2(
		rc + plh.z,
		rc * (1 - wgs84_epssq) + plh.z
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

	vec2 cisLambdaToTheM[nMax+1];
	cisLambdaToTheM[0] = vec2(1., 0.);
	cisLambdaToTheM[1] = cisLambda;

	for (int m=2; m <= nMax; ++m) {
		cisLambdaToTheM[m] = cplxmul(cisLambdaToTheM[m-1], cisLambda);
	}
<?
	
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
	
	float P[numTerms];	//Legendre function
	float dP[numTerms];	//Legendre function derivative
	P[0] = 1.;
	dP[0] = 0.;
<? 
	--	 First,	Compute the Gauss-normalized associated Legendre functions
	
	for n=1,nMax do 
		for m=0,n do
			local index = n * (n + 1) / 2 + m
			if n == m then
				local index1 = (n - 1) * n / 2 + m - 1
?>	P[<?=index?>] = cisPhiSph.x * P[<?=index1?>];
	dP[<?=index?>] = cisPhiSph.x * dP[<?=index1?>] + cisPhiSph.y * P[<?=index1?>];
<?			
			elseif n == 1 and m == 0 then
				local index1 = (n - 1) * n / 2 + m
?>	P[<?=index?>] = cisPhiSph.y * P[<?=index1?>];
	dP[<?=index?>] = cisPhiSph.y * dP[<?=index1?>] - cisPhiSph.x * P[<?=index1?>];
<?			
			elseif n > 1 and n ~= m then
				local index1 = (n - 2) * (n - 1) / 2 + m
				local index2 = (n - 1) * n / 2 + m
				if m > n - 2 then
?>	P[<?=index?>] = cisPhiSph.y * P[<?=index2?>];
	dP[<?=index?>] = cisPhiSph.y * dP[<?=index2?>] - cisPhiSph.x * P[<?=index2?>];
<?				
				else
					local k = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3))
?>	P[<?=index?>] = cisPhiSph.y * P[<?=index2?>] - <?=k?> * P[<?=index1?>];
	dP[<?=index?>] = cisPhiSph.y * dP[<?=index2?>] - cisPhiSph.x * P[<?=index2?>] - <?=k?> * dP[<?=index1?>];
<?
				end
			end
		end
	end


-- bake this into the B sum:
	
	-- Converts the Gauss-normalized associated Legendre
	-- functions to the Schmidt quasi-normalized version using pre-computed
	-- relation stored in the variable schmidtQuasiNorm
	
--	for n=1,nMax do
--		for m=0,n do
--			local index = n * (n + 1) / 2 + m
--	P[index] *= clnumber(schmidtQuasiNorm[index]);
--	dP[index] *= clnumber(-schmidtQuasiNorm[index]);
--		
--			-- The sign is changed since the new WMM routines use derivative with respect to latitude
--			-- insted of co-latitude
--		end
--	end
	
	-- end MAG_PcupLow
	-- end MAG_AssociatedLegendreFunction
	-- begin MAG_Summation 
?>

	float earthRadOverR = wgs84_re * invR;

	vec3 B = vec3(0., 0., 0.);
	{
		float earthRadOverRToTheN = earthRadOverR * earthRadOverR;
	
<? 
	for n=1,nMax do 
?>		earthRadOverRToTheN *= earthRadOverR;
<? 
		--.g .h .gt .ht == .xyzw
		-- then again, looks like I'm not doing any gt/ht calculations... 
		-- that means my reading is strictly 2020, right?
		
		--		    nMax  	(n+2) 	  n     m            m           m
		--		Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
		--						n=1      	      m=0   n            n           n 
		-- Equation 12 in the WMM Technical report.  Derivative with respect to radius.
?>		B.z -= 	earthRadOverRToTheN * (0.	
<?		for m=0,n do
			local index = (n * (n + 1) / 2 + m) 
?>			+ P[<?=index?>] * dot(cisLambdaToTheM[<?=m?>], vec2(<?=
				clnumber(wmm[n][m].g * (n + 1) * schmidtQuasiNorm[index])
			?>, <?=
				clnumber(wmm[n][m].h * (n + 1) * schmidtQuasiNorm[index])
			?>))
<?		end
?>		);
<?
		--		  1 nMax  (n+2)    n     m            m           m
		--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
		--				   n=1             m=0   n            n           n 
		-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.
?>		B.y += 	earthRadOverRToTheN * (0.	
<?		for m=1,n do
			local index = (n * (n + 1) / 2 + m) 
?>			+ P[<?=index?>] * dot(cisLambdaToTheM[<?=m?>], vec2(<?=
				clnumber(-wmm[n][m].h * m * schmidtQuasiNorm[index])
			?>, <?=
				clnumber(wmm[n][m].g * m * schmidtQuasiNorm[index])
			?>))
<?		end
?>		);
<?
		--		   nMax  (n+2) n     m            m           m
		--		Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
		--				   n=1         m=0   n            n           n 
		-- Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius.
?>		B.x -= 	earthRadOverRToTheN * (0.	
<?		for m=0,n do
			local index = (n * (n + 1) / 2 + m) 
?>			+ dP[<?=index?>] * dot(cisLambdaToTheM[<?=m?>], vec2(<?=
				clnumber(wmm[n][m].g * -schmidtQuasiNorm[index])
			?>, <?=
				clnumber(wmm[n][m].h * -schmidtQuasiNorm[index])
			?>))
<?		end 
?>		);
<?	end 
?>	}

	if (cisPhiSph.x < -1e-10 || cisPhiSph.x > 1e-10) {
		B.y /= cisPhiSph.x;
	} else {
		// Special calculation for component - By - at Geographic poles.
		// If the user wants to avoid using this function, please make sure that
		// the latitude is not exactly +/-90. An option is to make use the function
		// MAG_CheckGeographicPoles.
		// begin MAG_SummationSpecial	

		float PS[numTerms];
		PS[0] = 1.;

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
?>		PS[<?=n?>] = PS[<?=n-1?>];
<? 			else 
				local k = (((n - 1) * (n - 1)) - 1) / ((2 * n - 1) * (2 * n - 3))
?>		PS[<?=n?>] = cisPhiSph.y * PS[<?=n-1?>] - <?=clnumber(k)?> * PS[<?=n-2?>];
<? 
			end 

			--		  1 nMax  (n+2)    n     m            m           m
			--		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			--				   n=1             m=0   n            n           n 
			-- Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius.

?>		B.y += earthRadOverRToTheN * dot(cisLambda, vec2(<?=
			clnumber(-wmm[n][m].h * schmidtQuasiNorm3)
		?>, <?=
			clnumber(wmm[n][m].g * schmidtQuasiNorm3)
		?>)) * PS[<?=n?>];
<? 
		end 
		
		-- end MAG_SummationSpecial	
?>	}
	
	// end MAG_Summation 
	// end MAG_Geomag

	return B;
}

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

//for nMax == 12, numTerms == 91
const vec2 wmm[numTerms] = vec2[numTerms](
	vec2(0., 0.)	//n==0 is not used, but here for storage
<? 
for n=1,nMax do
	for m=0,n do
		local wmm_n_m = wmm[n][m]
?>	,vec2(<?=clnumber(wmm_n_m.g)?>, <?=clnumber(wmm_n_m.h)?>)
<?
	end
end 
?>);

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
	float r = length(xzp);
	float invR = 1. / r;
	vec2 cisPhiSph;
	cisPhiSph.y = xzp.y * invR;	// geocentric latitude sin & cos
	cisPhiSph.x = sqrt(1. - cisPhiSph.y * cisPhiSph.y);
	// longitude is the same 
	// end MAG_GeodeticToSpherical

	// begin MAG_Geomag
	// begin MAG_ComputeSphericalHarmonicVariables

	vec2 cisLambda = vec2(cos(plh.y), sin(plh.y));

	float earthRadOverR = wgs84_re * invR;

	vec2 cisLambdaToTheM[nMax+1];
	cisLambdaToTheM[0] = vec2(1., 0.);
	cisLambdaToTheM[1] = cisLambda;

	// looks like exp(i lambda)
	for (int m=2; m <= nMax; ++m) {
		cisLambdaToTheM[m] = cplxmul(cisLambdaToTheM[m-1], cisLambda);
	}
	
	// end MAG_ComputeSphericalHarmonicVariables
	// begin MAG_AssociatedLegendreFunction

	// begin MAG_PcupLow
	
	float x = cisPhiSph.y;
	float Pcup[numTerms];
	Pcup[0] = 1.;
	float dPcup[numTerms];
	dPcup[0] = 0.;

	// sin (geocentric latitude) - cisPhiSph.y
	float z = sqrt((1. - x) * (1. + x));

	//	 First,	Compute the Gauss-normalized associated Legendre functions
	for (int n=1; n <= nMax; ++n) {
		for (int m=0; m <= n; ++m) {
			int index = n * (n + 1) / 2 + m;
			if (n == m) {
				int index1 = (n - 1) * n / 2 + m - 1;
				Pcup[index] = z * Pcup[index1];
				dPcup[index] = z * dPcup[index1] + x * Pcup[index1];
			} else if (n == 1 && m == 0) { 
				int index1 = (n - 1) * n / 2 + m;
				Pcup[index] = x * Pcup[index1];
				dPcup[index] = x * dPcup[index1] - z * Pcup[index1];
			} else if (n > 1 && n != m) {
				int index1 = (n - 2) * (n - 1) / 2 + m;
				int index2 = (n - 1) * n / 2 + m;
				if (m > n - 2) {
					Pcup[index] = x * Pcup[index2];
					dPcup[index] = x * dPcup[index2] - z * Pcup[index2];
				} else {
					float k = float(((n - 1) * (n - 1)) - (m * m)) / float((2 * n - 1) * (2 * n - 3));
					Pcup[index] = x * Pcup[index2] - k * Pcup[index1];
					dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];
				}
			}
		}
	}
	// Compute the ration between the the Schmidt quasi-normalized associated Legendre
	// functions and the Gauss-normalized version. */

	float schmidtQuasiNorm[numTerms];
	schmidtQuasiNorm[0] = 1.;
	for (int n=1; n <= nMax; ++n) {
		int index = (n * (n + 1) / 2);
		int index1 = (n - 1) * n / 2;
		// for m = 0
		schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (2 * n - 1) / float(n);

		for (int m=1; m <= n; ++m) {
			int index = (n * (n + 1) / 2 + m);
			int index1 = (n * (n + 1) / 2 + m - 1);
			schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * sqrt( (float(n - m + 1) * (m == 1 ? 2. : 1.)) / float(n + m));
		}
	}

	// Converts the Gauss-normalized associated Legendre
	// functions to the Schmidt quasi-normalized version using pre-computed
	// relation stored in the variable schmidtQuasiNorm */

	for (int n=1; n <= nMax; ++n) {
		for (int m=0; m <= n; ++m) {
			int index = (n * (n + 1) / 2 + m);
			Pcup[index] = Pcup[index] * schmidtQuasiNorm[index];
			dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index];
			// The sign is changed since the new WMM routines use derivative with respect to latitude
			// insted of co-latitude */
		}
	}
	
	// end MAG_PcupLow
	// end MAG_AssociatedLegendreFunction
	// begin MAG_Summation 

	vec3 B = vec3(0., 0., 0.);
	{
		float earthRadOverRToTheN = earthRadOverR * earthRadOverR;
		for (int n=1; n <= nMax; ++n) {
			earthRadOverRToTheN *= earthRadOverR;
			for (int m=0; m <= n; ++m) {
				int index = (n * (n + 1) / 2 + m);
				
				//.g .h .gt .ht == .xyzw
				// then again, looks like I'm not doing any gt/ht calculations... 
				// that means my reading is strictly 2020, right?
				vec2 wmm_n_m = wmm[index];	//[n][m]

				//		    nMax  	(n+2) 	  n     m            m           m
				//		Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
				//						n=1      	      m=0   n            n           n  */
				// Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
				B.z -=
					earthRadOverRToTheN
					* (
						wmm_n_m.x * cisLambdaToTheM[m].x
						+ wmm_n_m.y * cisLambdaToTheM[m].y
					)
					* (n + 1) * Pcup[index];

				//		  1 nMax  (n+2)    n     m            m           m
				//		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
				//				   n=1             m=0   n            n           n  */
				// Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
				B.y += 
					earthRadOverRToTheN
					* (
						wmm_n_m.x * cisLambdaToTheM[m].y
						- wmm_n_m.y * cisLambdaToTheM[m].x
					)
					* float(m) * Pcup[index];
				//		   nMax  (n+2) n     m            m           m
				//		Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
				//				   n=1         m=0   n            n           n  */
				// Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */
				B.x -=
					earthRadOverRToTheN *
					(
						wmm_n_m.x * cisLambdaToTheM[m].x
						+ wmm_n_m.y * cisLambdaToTheM[m].y
					)
					* dPcup[index];
			}
		}
	}

	if (cisPhiSph.x < -1e-10 || cisPhiSph.x > 1e-10) {
		B.y /= cisPhiSph.x;
	} else {
		// Special calculation for component - By - at Geographic poles.
		// If the user wants to avoid using this function, please make sure that
		// the latitude is not exactly +/-90. An option is to make use the function
		// MAG_CheckGeographicPoles.
		// begin MAG_SummationSpecial	

		float PcupS[numTerms];
		PcupS[0] = 1.;
		float schmidtQuasiNorm1 = 1.;

		B.y = 0.;

		float earthRadOverRToTheN = earthRadOverR * earthRadOverR;
		for (int n=1; n <= nMax; ++n) {
			earthRadOverRToTheN *= earthRadOverR;
			
			//Compute the ration between the Gauss-normalized associated Legendre
			// functions and the Schmidt quasi-normalized version. This is equivalent to
			// sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)! */
			const int m = 1;
			int index = n * (n + 1) / 2 + m;
			vec2 wmm_n_m = wmm[index];
			
			float schmidtQuasiNorm2 = schmidtQuasiNorm1 * float(2 * n - 1) / float(n);
			float schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt( float(n * 2) / float(n + 1));
			float schmidtQuasiNorm1 = schmidtQuasiNorm2;
			if (n == 1) {
				PcupS[n] = PcupS[n-1];
			} else {
				float k = float(((n - 1) * (n - 1)) - 1) / float((2 * n - 1) * (2 * n - 3));
				PcupS[n] = cisPhiSph.y * PcupS[n - 1] - k * PcupS[n - 2];
			}

			//		  1 nMax  (n+2)    n     m            m           m
			//		By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
			//				   n=1             m=0   n            n           n  */
			// Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
			B.y += 
				earthRadOverRToTheN
				* (
					wmm_n_m.x * cisLambdaToTheM[1].y
					- wmm_n_m.y * cisLambdaToTheM[1].x
				) * PcupS[n] * schmidtQuasiNorm3;
		}

		// end MAG_SummationSpecial	
	}
	
	// end MAG_Summation 
	// end MAG_Geomag

	return B;
}

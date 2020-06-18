// This program implements conditional or exact ML estimation
// of a mixed-frequency VAR model without a constant term
// by applying the EM algorithm and a quasi-Newton method (BFGS) sequentially.
// The user must demean the input series before estimation.

#include <oxstd.h>
#include <oxfloat.h>
#include <packages/ssfpack/ssfpack.h>
#import <maximize>

class VAREM
{
	GetSSF(const iP, const cN_1, const vY, const mPhi, const mSigma,
		const amA, const amB, const amC, const amD);
	prediction(const vS_u, const mP_u, const mA, const mB, const avS_p, const amP_p);
	updating(const vS_p, const mP_p, const vE, const mC, const mD, const avS_u, const amP_u);
	KF(const mY, const mPhi, const mSigma, const amS_p, const amVechP_p, const avLnL);
	smoothing(const vS_p, const mP_p, const vE, const vR, const mR, const mA, const mB, const mC, const mD,
		const avS_s, const amP_s, const avR, const amR);
	KS(const mY, const mPhi, const mSigma, const mS_p, const mVechP_p,
		const amS_s, const amVechP_s, const amVecP1_s);
	EM(const mY, const mPhi, const mSigma, const amPhi, const amSigma);
	loglikelihood(const mY, const mPhi, const mSigma);
	estimate(const mY, const iP, const amPhi, const amSigma, const adLnL);
}

VAREM::GetSSF(const iP, const cN_1, const vY, const mPhi, const mSigma,
	const amA, const amB, const amC, const amD)
{
	// Transition equation

	decl cn = rows(vY);
	decl ma;
	if (iP <= 5)
	{
		ma = (mPhi ~ zeros(cn, (5 - iP) * cn)) | (unit(4 * cn) ~ zeros(4 * cn, cn));
	}
	else
	{
		ma = mPhi | (unit((iP - 1) * cn) ~ zeros((iP - 1) * cn, cn));
	}
	amA[0] = ma;
	decl cm = rows(ma);
	amB[0] = choleski(mSigma) | zeros(cm - cn, cn);
	

	// Measurement equation

	decl cn_2 = cn - cN_1;
	decl mc_0 = ((1 / 3) * unit(cN_1) ~ zeros(cN_1, cn_2)) | (zeros(cn_2, cN_1) ~ unit(cn_2));
	decl mc_1 = ((2 / 3) * unit(cN_1) ~ zeros(cN_1, cn_2)) | zeros(cn_2, cn);
	decl mc_2 = (          unit(cN_1) ~ zeros(cN_1, cn_2)) | zeros(cn_2, cn);
	decl mc_3 = ((2 / 3) * unit(cN_1) ~ zeros(cN_1, cn_2)) | zeros(cn_2, cn);
	decl mc_4 = ((1 / 3) * unit(cN_1) ~ zeros(cN_1, cn_2)) | zeros(cn_2, cn);
	decl mc = mc_0 ~ mc_1 ~ mc_2 ~ mc_3 ~ mc_4;
	if (iP > 5)
	{
		mc = mc ~ zeros(cn, (iP - 5) * cn);
	}
	decl md = zeros(cn, cn);
	if (isnan(vY))
	{
		mc[:cN_1 - 1][] = 0;
		md[:cN_1 - 1][] = unit(cN_1);
	}
	amC[0] = mc;
	amD[0] = md;
}

VAREM::prediction(const vS_u, const mP_u, const mA, const mB, const avS_p, const amP_p)
{
	avS_p[0] = mA * vS_u;
	amP_p[0] = mA * mP_u * mA' + mB * mB';
}

VAREM::updating(const vS_p, const mP_p, const vE, const mC, const mD, const avS_u, const amP_u)
{
	decl cm = rows(vS_p);
	decl mk = mP_p * mC' * invertsym(mC * mP_p * mC' + mD * mD');
	avS_u[0] = vS_p + mk * vE;
	amP_u[0] = (unit(cm) - mk * mC) * mP_p;
}

VAREM::KF(const mY, const mPhi, const mSigma, const amS_p, const amVechP_p, const avLnL)
{
	decl cn = rows(mY);
	decl cn_1 = rows(selectr(mY));
	decl ct = columns(mY);
	decl ip = columns(mPhi) / cn;
	decl cm = max(5, ip) * cn;
	decl ms_p = zeros(cm, ct);
	decl mvechP_p = zeros(cm * (cm + 1)/ 2, ct);
	decl vlnL = zeros(ct, 1);
	decl vs_u = zeros(cm, 1);
	decl mp_u = zeros(cm, cm);
	decl t;
	for (t = 0; t < ct; ++t)
	{
		decl vy = mY[][t];
		decl ma, mb, mc, md;
		GetSSF(ip, cn_1, vy, mPhi, mSigma, &ma, &mb, &mc, &md);

		decl vs_p, mp_p;
		prediction(vs_u, mp_u, ma, mb, &vs_p, &mp_p);
		ms_p[][t] = vs_p;
		mvechP_p[][t] = vech(mp_p);

		if (isnan(vy))
		{
			vy[:cn_1 - 1] = zeros(cn_1, 1);
		}
		decl ve = vy - mc * vs_p;
		decl msigma_ee = mc * mp_p * mc' + md * md';
		decl dlnL = - (cn / 2) * log(2 * M_PI)
			- (1 / 2) * log(determinant(msigma_ee))
			- (1 / 2) * ve' * invertsym(msigma_ee) * ve;
		vlnL[t] = dlnL;

		updating(vs_p, mp_p, ve, mc, md, &vs_u, &mp_u);
	}
	amS_p[0] = ms_p;
	amVechP_p[0] = mvechP_p;
	avLnL[0] = vlnL;
}

VAREM::smoothing(const vS_p, const mP_p, const vE, const vR, const mR, const mA, const mB, const mC, const mD,
	const avS_s, const amP_s, const avR, const amR)
{
	decl cm = rows(vS_p);
	decl mk = mP_p * mC' * invertsym(mC * mP_p * mC' + mD * mD');
	decl ml = mA * (unit(cm) - mk * mC);
	decl vr = mC' * invertsym(mC * mP_p * mC' + mD * mD') * vE + ml' * vR;
	decl mr = mC' * invertsym(mC * mP_p * mC' + mD * mD') * mC + ml' * mR * ml;
	avS_s[0] = vS_p + mP_p * vr;
	amP_s[0] = mP_p - mP_p * mr * mP_p;
	avR[0] = vr;
	amR[0] = mr;
}

VAREM::KS(const mY, const mPhi, const mSigma, const mS_p, const mVechP_p,
	const amS_s, const amVechP_s, const amVecP1_s)
{
	decl cn = rows(mY);
	decl cn_1 = rows(selectr(mY));
	decl ct = columns(mY);
	decl ip = columns(mPhi) / cn;
	decl cm = max(5, ip) * cn;
	decl ms_s = zeros(cm, ct);
	decl mvechP_s = zeros(cm * (cm + 1)/ 2, ct);
	decl mvecP1_s = zeros(cm * cm, ct);
	decl vr = zeros(cm, 1);
	decl mr = zeros(cm, cm);
	decl t;
	for (t = ct - 1; t >= 0; --t)
	{
		decl vy = mY[][t];
		decl ma, mb, mc, md;
		GetSSF(ip, cn_1, vy, mPhi, mSigma, &ma, &mb, &mc, &md);

		decl vs_p = mS_p[][t];
		decl mp_p = unvech(mVechP_p[][t]);
		if (isnan(vy))
		{
			vy[:cn_1 - 1] = zeros(cn_1, 1);
		}
		decl ve = vy - mc * vs_p;
		decl vs_s, mp_s;
		smoothing(vs_p, mp_p, ve, vr, mr, ma, mb, mc, md, &vs_s, &mp_s, &vr, &mr);
		ms_s[][t] = vs_s;
		mvechP_s[][t] = vech(mp_s);
		if (t > 0)
		{
			decl mlp_p = unvech(mVechP_p[][t - 1]);
			decl vly = mY[][t - 1];
			decl mla, mlb, mlc, mld;
			GetSSF(ip, cn_1, vly, mPhi, mSigma, &mla, &mlb, &mlc, &mld);

			decl mlk = mlp_p * mlc' * invertsym(mlc * mlp_p * mlc' + mld * mld');
			decl mll = mla * (unit(cm) - mlk * mlc);
			decl mp1_s = (unit(cm) - mp_p * mr) * mll * mlp_p;
			mvecP1_s[][t] = vec(mp1_s);
		}
	}
	amS_s[0] = ms_s;
	amVechP_s[0] = mvechP_s;
	amVecP1_s[0] = mvecP1_s;
}

VAREM::EM(const mY, const mPhi, const mSigma, const amPhi, const amSigma)
{
	// E step

	decl ms_p, mvechP_p, vlnL;
	KF(mY, mPhi, mSigma, &ms_p, &mvechP_p, &vlnL);
	decl ms_s, mvechP_s, mvecP1_s;
	KS(mY, mPhi, mSigma, ms_p, mvechP_p, &ms_s, &mvechP_s, &mvecP1_s);
	
	// M step

	decl cn = rows(mY);
	decl ct = columns(mY);
	decl ip = columns(mPhi) / cn;
	decl cm = max(5, ip) * cn;
	decl mls_s = lag0(ms_s', 1)';
	decl mp = unvech(sumr(mvechP_s) / ct);
	decl mp1 = shape(sumr(mvecP1_s) / ct, cm, cm);
	decl mlp = unvech(sumr(lag0(mvechP_s', 1)') / ct);
	decl mm_00 = mp + ms_s * ms_s' / ct;
	decl mm_01 = mp1 + ms_s * mls_s' / ct;
	decl mm_10 = mm_01';
	decl mm_11 = mlp + mls_s * mls_s' / ct;
	decl mf = unit(cn) ~ zeros(cn, (ip - 1) * cn);
	decl mg = unit(ip * cn);
	if (ip <= 5)
	{
		mf = unit(cn) ~ zeros(cn, 4 * cn);
		mg = unit(ip * cn) ~ zeros(ip * cn, (5 - ip) * cn);
	}
	decl mphi = mf * mm_01 * mg' * invertsym(mg * mm_11 * mg');
	decl msigma = mf * (mm_00 - mm_01 * mg' * invertsym(mg * mm_11 * mg') * mg * mm_10) * mf';
	amPhi[0] = mphi;
	amSigma[0] = msigma;
}

VAREM::loglikelihood(const mY, const mPhi, const mSigma)
{
	decl ms_p, mvechP_p, vlnL;
	KF(mY, mPhi, mSigma, &ms_p, &mvechP_p, &vlnL);
	return meanc(vlnL);
}

VAREM::estimate(const mY, const iP, const amPhi, const amSigma, const adLnL)
{
	decl cn = rows(mY);
	decl mphi_old = zeros(cn, iP * cn);
	decl msigma_old = unit(cn);
	decl dlnL_old = loglikelihood(mY, mphi_old, msigma_old);
	decl mphi, msigma, dlnL, d;
	do
	{
		EM(mY, mphi_old, msigma_old, &mphi, &msigma);
		dlnL = loglikelihood(mY, mphi, msigma);
		d = dlnL - dlnL_old;
		mphi_old = mphi;
		msigma_old = msigma;
		dlnL_old = dlnL;
		print(dlnL ~ d);
//		print(mphi, msigma);
	} while (d > 10 ^ -8);  // Convergence criterion
	amPhi[0] = mphi;
	amSigma[0] = msigma;
	adLnL[0] = dlnL;
}

class VARQN
{
	decl m_mY, m_iP;
	VARQN(const mY, const iP);
	GetVAR(const vTheta, const cN, const iP, const amPi, const amSigma_w);
	GetSSF(const mPi, const mSigma_w, const amPhi, const amOmega, const amSigma);
	loglikelihood(const vTheta, const adLnL, const avScore, const amHess);
	estimate(const mY, const iP, const mPi_0, const mSigma_w0, const amPi, const amSigma_w, const adLnL);
	smoothing(const mY, const mPi, const mSigma_w);
}

VARQN::VARQN(const mY, const iP)
{
	m_mY = mY;
	m_iP = iP;
}

VARQN::GetVAR(const vTheta, const cN, const iP, const amPi, const amSigma_w)
{
	decl ck = iP * cN * cN + cN * (cN + 1) / 2;
	amPi[0] = shape(vTheta[:iP * cN * cN - 1], cN, iP * cN);
	decl mcholeski = lower(unvech(vTheta[iP * cN * cN: ck - 1]));
	amSigma_w[0] = mcholeski * mcholeski';
}

VARQN::GetSSF(const mPi, const mSigma_w, const amPhi, const amOmega, const amSigma)
{
	decl my = m_mY;
	decl cn = rows(my);
	decl cn_1 = rows(selectr(my));
	decl cn_2 = cn - cn_1;
	decl ip = columns(mPi) / cn;

	// Transition equation

	decl ma;
	if (ip <= 5)
	{
		ma = (mPi ~ zeros(cn, (5 - ip) * cn)) | (unit(4 * cn) ~ zeros(4 * cn, cn));
	}
	else
	{
		ma = mPi | (unit((ip - 1) * cn) ~ zeros((ip - 1) * cn, cn));
	}
	decl cm = rows(ma);
	decl mb = choleski(mSigma_w) | zeros(cm - cn, cn);

	// Measurement equation

	decl mc_0 = ((1 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2)) | (zeros(cn_2, cn_1) ~ unit(cn_2));
	decl mc_1 = ((2 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2)) | zeros(cn_2, cn);
	decl mc_2 = (          unit(cn_1) ~ zeros(cn_1, cn_2)) | zeros(cn_2, cn);
	decl mc_3 = ((2 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2)) | zeros(cn_2, cn);
	decl mc_4 = ((1 / 3) * unit(cn_1) ~ zeros(cn_1, cn_2)) | zeros(cn_2, cn);
	decl mc = mc_0 ~ mc_1 ~ mc_2 ~ mc_3 ~ mc_4;
	if (ip > 5)
	{
		mc = mc ~ zeros(cn, (ip - 5) * cn);
	}
	decl md = zeros(cn, cn);
	decl me = mb | md;
	amPhi[0] = ma | mc;
	amOmega[0] = me * me';
//	amSigma[0] = shape(invert(unit(cm ^ 2) - ma ** ma) * vec(mb * mb'), cm, cm) | zeros(1, cm);  // exact ML;
	amSigma[0] = mb * mb' | zeros(1, cm);  // conditional ML;
}

VARQN::loglikelihood(const vTheta, const adLnL, const avScore, const amHess)
{
	decl my = m_mY;
	decl cn = rows(my);
	decl ct = columns(my);
	decl ip = m_iP;
	decl mpi, msigma_w;
	GetVAR(vTheta, cn, ip, &mpi, &msigma_w);
	decl mphi, momega, msigma;
	GetSSF(mpi, msigma_w, &mphi, &momega, &msigma);
	decl dlnL, dvar;
	SsfLik(&dlnL, &dvar, my, mphi, momega, msigma);
	adLnL[0] = dlnL / ct;
	return 1;
}

VARQN::estimate(const mY, const iP, const mPi_0, const mSigma_w0, const amPi, const amSigma_w, const adLnL)
{
	decl vtheta = vec(mPi_0) | vech(choleski(mSigma_w0));
	decl dfunc;
	MaxControl(-1, 1);
	MaxBFGS(loglikelihood, &vtheta, &dfunc, 0, 1);
	decl cn = rows(mY);
	decl ct = columns(mY);
//	decl mHess;
//	Num2Derivative(loglikelihood, vtheta, &mHess);
//	print(sqrt(diagonal(invertsym(-ct * mHess)))');
	decl mpi, msigma_w;
	GetVAR(vtheta, cn, iP, &mpi, &msigma_w);
	amPi[0] = mpi;
	amSigma_w[0] = msigma_w;
	adLnL[0] = dfunc;
}

VARQN::smoothing(const mY, const mPi, const mSigma_w)
{
	decl mphi, momega, msigma;
	GetSSF(mPi, mSigma_w, &mphi, &momega, &msigma);
	decl ms_s;
	decl mSmo = SsfMomentEst(ST_SMO, &ms_s, mY, mphi, momega, msigma);
	return ms_s;
}

main()
{
	decl mY = loadmat("bciq1m4.xls")';
	decl cn = rows(mY);
	decl ct = columns(mY);
	decl ip = 1;
//	for (ip = 1; ip <= 12; ++ip)  // Model selection
	{
		decl mphi = zeros(cn, ip * cn);
		decl msigma = unit(cn);
		decl dlnL;

		decl varemobj = new VAREM();
		varemobj.estimate(mY, ip, &mphi, &msigma, &dlnL);  // Comment out if EM is unnecessary.
		delete varemobj;

		decl varqnobj = new VARQN(mY, ip);
		varqnobj.estimate(mY, ip, mphi, msigma, &mphi, &msigma, &dlnL);  // Comment out if QN is unnecessary.
		decl ms_s = varqnobj.smoothing(mY, mphi, msigma);
		delete varqnobj;

		dlnL = ct * dlnL;
		decl ck = ip * cn * cn;
		decl dAIC = (dlnL - ck) / ct;
		decl dBIC = (dlnL - ck * log(ct) / 2) / ct;
		print(mphi, msigma);
		print(ct);
		print(ip ~ dlnL ~ dAIC ~ dBIC);

		decl vindex = ms_s[0][]';
		savemat("c:/Program Files/OxMetrics4/Ox/myprg_ox/04a/temp.xls", vindex);
//		print(vindex);
	}
}
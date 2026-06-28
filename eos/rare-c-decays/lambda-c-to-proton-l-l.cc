/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Ahmet Kokulu
 * Copyright (c) 2019,2021 Danny van Dyk
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2026 Dominik Suelmann
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/form-factors/baryonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-c-decays/lambda-c-to-proton-l-l.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using namespace std::literals::string_literals;
    using std::norm;

    namespace lambda_c_to_proton_l_l
    {
        struct Amplitudes
        {
                double aU11;
                double aL11;
                double aU22;
                double aL22;
                double aP12;
                double aS22;

                double alpha;
                double beta_l;
                double beta_l_squared;
                double four_m_l_squared;
                double q2;
                double norm_squared;
        };

        struct AngularObservables
        {
                std::array<double, 3> _k;

                AngularObservables(const Amplitudes & a)
                {
                    using std::conj;
                    using std::imag;
                    using std::norm;
                    using std::real;
                    using std::sqrt;

                    double beta_l           = a.beta_l;
                    double beta_l_squared   = a.beta_l_squared;
                    double four_m_l_squared = a.four_m_l_squared;
                    double q2               = a.q2;
                    double aU11U22          = a.aU11 + a.aU22;
                    double aULS             = a.aU11 + a.aL11 + a.aS22;

                    // cf. [GHM:2021A], eq. (20), p. 6.

                    // K_{1ss}
                    _k[0] = a.norm_squared * (q2 * beta_l_squared * (aU11U22 / 2.0 + a.aL11 + a.aL22) + four_m_l_squared * aULS);

                    // K_{1cc}
                    _k[1] = a.norm_squared * (q2 * beta_l_squared * aU11U22 + four_m_l_squared * aULS);

                    // K_{1c}
                    _k[2] = a.norm_squared * (-2.0 * q2 * beta_l * a.aP12);
                }

                AngularObservables(const std::array<double, 3> & k) :
                    _k(k)
                {
                }

                inline double
                k1ss() const
                {
                    return _k[0];
                }

                inline double
                k1cc() const
                {
                    return _k[1];
                }

                inline double
                k1c() const
                {
                    return _k[2];
                }

                inline double
                decay_width() const
                {
                    return 2.0 * k1ss() + k1cc();
                }

                inline double
                a_fb_leptonic() const
                {
                    return 3.0 / 2.0 * k1c() / decay_width();
                }

                inline double
                a_fb_leptonic_num() const
                {
                    return 3.0 / 2.0 * k1c();
                }

                inline double
                f_l() const
                {
                    return (2.0 * k1ss() - k1cc()) / decay_width();
                }

                inline double
                f_l_num() const
                {
                    return (2.0 * k1ss() - k1cc());
                }

                inline double
                d2gamma(const double & c_lep) const
                {
                    const double c2_lep = c_lep * c_lep;
                    const double s2_lep = 1.0 - c2_lep;
                    const double s_lep  = std::sqrt(s2_lep);

                    // cf. [GHM:2021A], p. 6, eq. (9)
                    return 3.0 / 2.0 * (k1ss() * s2_lep + k1cc() * c2_lep + k1c() * c_lep);
                }
        };
    } // namespace lambda_c_to_proton_l_l

    template <> struct Implementation<LambdaCToProtonLeptonLepton>
    {
            std::shared_ptr<Model> model;
            std::shared_ptr<Model> model_SM;

            Parameters parameters;

            UsedParameter hbar;
            UsedParameter tau_Lambda_c;

            UsedParameter g_fermi;
            UsedParameter alpha_e;

            LeptonFlavorOption opt_l;
            UsedParameter      m_l;

            UsedParameter m_Lambda_c;
            UsedParameter m_proton;
            UsedParameter m_rho;
            UsedParameter m_omega;
            UsedParameter m_phi;
            UsedParameter gamma_rho;
            UsedParameter gamma_omega;
            UsedParameter gamma_phi;

            // resonance parameters
            UsedParameter a_rho;
            UsedParameter a_omega;
            UsedParameter a_phi;
            UsedParameter delta_rho;
            UsedParameter delta_omega_m_rho;
            UsedParameter delta_phi_m_rho;

            BooleanOption opt_cp_conjugate;

            UsedParameter mu;
            UsedParameter ms;
            UsedParameter mc;

            // CKM parameters
            UsedParameter abs_Vud;
            UsedParameter abs_Vus;
            UsedParameter abs_Vcd;
            UsedParameter abs_Vcs;
            UsedParameter arg_Vcd;
            UsedParameter arg_Vcs;

            std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

            static const std::vector<OptionSpecification> options;

            Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
                model(Model::make(o.get("model"_ok, "WET"), p, o)),
                model_SM(Model::make("SM", p, o)),
                parameters(p),
                hbar(p["QM::hbar"], u),
                tau_Lambda_c(p["life_time::Lambda_c"], u),
                g_fermi(p["WET::G_Fermi"], u),
                alpha_e(p["QED::alpha_e(m_c)"], u),
                opt_l(o, options, "l"_ok),
                m_l(p["mass::" + opt_l.str()], u),
                m_Lambda_c(p["mass::Lambda_c"], u),
                m_proton(p["mass::proton"], u),
                m_rho(p["mass::rho^0"], u),
                m_omega(p["mass::omega"], u),
                m_phi(p["mass::phi"], u),
                gamma_rho(p["gamma::rho^0"], u),
                gamma_omega(p["gamma::omega"], u),
                gamma_phi(p["gamma::phi"], u),
                a_rho(p["Lambda_c->proton::res_a_rho"], u),
                a_omega(p["Lambda_c->proton::res_a_omega"], u),
                a_phi(p["Lambda_c->proton::res_a_phi"], u),
                delta_rho(p["Lambda_c->proton::res_delta_rho"], u),
                delta_omega_m_rho(p["Lambda_c->proton::res_delta_omega_m_rho"], u),
                delta_phi_m_rho(p["Lambda_c->proton::res_delta_phi_m_rho"], u),
                opt_cp_conjugate(o, options, "cp-conjugate"_ok),
                mu(p["uc::mu"], u),
                ms(p["mass::s(2GeV)"], u),
                mc(p["mass::c"], u),
                abs_Vud(p["CKM::abs(V_ud)"], u),
                abs_Vus(p["CKM::abs(V_us)"], u),
                abs_Vcd(p["CKM::abs(V_cd)"], u),
                abs_Vcs(p["CKM::abs(V_cs)"], u),
                arg_Vcd(p["CKM::arg(V_cd)"], u),
                arg_Vcs(p["CKM::arg(V_cs)"], u),
                form_factors(FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_c->neutron::" + o.get("form-factors"_ok, "BMRvD2022"), p, o))
            {
                u.uses(*form_factors);
                u.uses(*model);
            }

            const complex<double>
            normalization(const double & q2) const
            {
                double lam = lambda(m_Lambda_c * m_Lambda_c, m_proton * m_proton, q2);

                if ((lam <= 0) || (q2 <= 4.0 * m_l * m_l))
                {
                    return 0.0;
                }

                // alpha_e = 0.007501720229738859
                return g_fermi() * 0.007501720229738859 * std::sqrt(std::sqrt(1.0 - 4.0 * m_l * m_l / q2))
                       * std::sqrt(std::sqrt(lam) / 3.0 / power_of<11>(2) / power_of<5>(M_PI) / power_of<3>(m_Lambda_c()));
            }

            const complex<double>
            C9eff(const double & q2, const eos::WilsonCoefficients<eos::wc::UC> & c) const
            {
                complex<double> Lzero(5. / 3. + log(mu * mu / q2), M_PI);
                double          x_mc = 4 * mc * mc / q2;
                complex<double> fak_mc;
                if (x_mc < 1)
                {
                    fak_mc = log((1 + sqrt(1 - x_mc)) / (1 - sqrt(1 - x_mc))) - complex<double>(0., M_PI);
                }
                else
                {
                    fak_mc = 2 * std::atan(1 / sqrt(x_mc - 1));
                }
                complex<double> Lmc = 5. / 3. + log(mu * mu / mc) + x_mc - 0.5 * (2 + x_mc) * std::sqrt(std::abs(1 - x_mc)) * fak_mc;

                double          x_ms = 4 * ms * ms / q2;
                complex<double> fak_ms;
                if (x_ms < 1)
                {
                    fak_ms = log((1 + sqrt(1 - x_ms)) / (1 - sqrt(1 - x_ms))) - complex<double>(0., M_PI);
                }
                else
                {
                    fak_ms = 2 * std::atan(1 / sqrt(x_ms - 1));
                }
                complex<double> Lms = 5. / 3. + log(mu * mu / ms) + x_ms - 0.5 * (2 + x_ms) * std::sqrt(std::abs(1 - x_ms)) * fak_ms;

                complex<double> ret;
                ret = c.c9_i(2)
                      + c._alpha_s / (4 * M_PI)
                                * (8. / 27. * c.c1_i(1) + 2. / 9. * c.c2_i(1) - 8. / 9. * c.c3_i(1) - 32. / 27. * c.c4_i(1) - 128. / 9. * c.c5_i(1) - 512. / 27. * c.c6_i(1)
                                   + Lmc * (28. / 9. * c.c3_i(1) + 16. / 27. * c.c4_i(1) + 304. / 9. * c.c5_i(1) + 256. / 27. * c.c6_i(1))
                                   + Lms * (-4. / 3. * c.c3_i(1) - 40. / 3. * c.c5_i(1))
                                   + Lzero * (16. / 9. * c.c3_i(1) + 16. / 27. * c.c4_i(1) + 184. / 9. * c.c5_i(1) + 256. / 27. * c.c6_i(1)));

                complex<double> C8_eff = c.c8_i(1) + c.c3_i(0) - 1. / 6. * c.c4_i(0) + 20. * c.c5_i(0) - 10. / 3. * c.c6_i(0);

                double          Q2 = q2 / (mc * mc);
                complex<double> F_89 =
                        (-16. * M_PI * M_PI / 27. * (4. - Q2) / ((1. - Q2) * (1. - Q2) * (1. - Q2) * (1. - Q2)) + 16. / 9. * (5. - 2. * Q2) / ((1. - Q2) * (1. - Q2))
                         + 32. / 9. * std::sqrt(4. - Q2) / std::sqrt(Q2) / ((1. - Q2) * (1. - Q2) * (1. - Q2)) * (4. + 3. * Q2 - Q2 * Q2) * std::asin(std::sqrt(Q2) / 2.)
                         + 64. / 3. * (4. - Q2) / ((1. - Q2) * (1. - Q2) * (1. - Q2) * (1. - Q2)) * std::asin(std::sqrt(Q2) / 2.) * std::asin(std::sqrt(Q2) / 2.)
                         + 32. / 9. * 1. / (1. - Q2) * std::log(Q2));

                complex<double> q2_c(q2, 0.);
                complex<double> I(0., 1.);
                complex<double> F_92_d =
                        ((16.850454731592844 - 1.6541601433384439 * I) - (0.09239330957393718 - 0.9828824718848296 * I) / (q2 * q2)
                         - (3.2505185476820753 + 9.517218005264239 * I) / q2 - (3.8840518779093536 - 1.915097525241604 * I) * q2
                         + (0.15187990867462747 - 0.34719611566798503 * I) * q2 * q2 + (3.2928401123637885 + 4.341576642611662 * I) * std::log(mu)
                         + ((0.0004041743729919311 - 0.00018938373677914851 * I) * std::log(mu)) / (q2 * q2)
                         - ((0.02906629754068944 - 0.020602741706600855 * I) * std::log(mu)) / q2 - (0.014111270976663197 - 0.013977145126954837 * I) * q2 * std::log(mu)
                         - 0.002289471791020287 * q2 * q2 * std::log(mu) + 0.4938261662214249 * std::log(mu) * std::log(mu)
                         - (7.977556163485848 - 1.3865595256182361 * I) * std::log(-q2_c) - ((0.0020702052990949613 - 0.13523901011279155 * I) * std::log(-q2_c)) / (q2 * q2)
                         + ((4.441100060097692 - 1.8216265342986375 * I) * std::log(-q2_c)) / q2 + (1.825860818894409 + 0.20764261026655212 * I) * q2 * std::log(-q2_c)
                         + (0.015315601748722353 + 0.002645487297212485 * I) * std::log(mu) * std::log(-q2_c)
                         - ((0. + 0.00013966961204803138 * I) * std::log(mu) * std::log(-q2_c)) / (q2 * q2)
                         - ((0.007779866128356353 - 0.0012317190034184576 * I) * std::log(mu) * std::log(-q2_c)) / q2
                         - (0.004448227613770954 + 0.0014304210693156473 * I) * q2 * std::log(mu) * std::log(-q2_c)
                         + (0.2654091420687322 - 2.417837791888905 * I) * std::log(-q2_c) * std::log(-q2_c)
                         - ((0.0030052343022420034 - 0.010123939239126046 * I) * std::log(-q2_c) * std::log(-q2_c)) / (q2 * q2)
                         + ((0.7539928404989259 + 0.12771928827451795 * I) * std::log(-q2_c) * std::log(-q2_c)) / q2
                         - (0.0014595576227236183 - 0.00463707219667801 * I) * std::log(mu) * std::log(-q2_c) * std::log(-q2_c)
                         - ((0.0010842440443838643 + 0.0003891352408633328 * I) * std::log(mu) * std::log(-q2_c) * std::log(-q2_c)) / q2
                         - (11.524567190326126 + 7.729142929989321 * I) * std::log(q2) - ((0.07627975814470198 - 0.13207926572744955 * I) * std::log(q2)) / (q2 * q2)
                         - ((5.26716180778514 + 2.2750628876413237 * I) * std::log(q2)) / q2 - (1.0673195590787947 + 1.2242703742442957 * I) * q2 * std::log(q2)
                         - (1.3965409742749544 - 0.00652508856088367 * I) * std::log(mu) * std::log(q2)
                         + ((0. + 0.0001094532346470768 * I) * std::log(mu) * std::log(q2)) / (q2 * q2)
                         - ((0.005337620318113501 - 0.005581477052076378 * I) * std::log(mu) * std::log(q2)) / q2
                         - (0.004449309854525112 - 0.0014301039268257902 * I) * q2 * std::log(mu) * std::log(q2)
                         + (2.7606095016443657 - 1.5519565733309244 * I) * std::log(q2) * std::log(q2)
                         - ((0.008066520393461691 - 0.008746618593736705 * I) * std::log(q2) * std::log(q2)) / (q2 * q2)
                         - ((0.5850242438810518 + 0.7224819982868034 * I) * std::log(q2) * std::log(q2)) / q2
                         - (0.0014591126757388096 + 0.004638490089092232 * I) * std::log(mu) * std::log(q2) * std::log(q2)
                         - ((0.0010846907950053178 - 0.0003890119975143421 * I) * std::log(mu) * std::log(q2) * std::log(q2)) / q2);
                complex<double> F_92_s =
                        ((131.3105122735358 + 32.21025602851325 * I) - (7.571659984976852 - 3.7429791914883506 * I) / (-0.06760000000000001 + q2)
                         - (0.21927851191580947 + 57.140118516713265 * I) * q2 - ((80.18453993786741 - 54.40263417725573 * I) * q2) / (-0.06760000000000001 + q2)
                         + (2.5875672915583285 - 1.18252495315284 * I) * q2 * q2 - ((29.646332913813794 + 35.49491507299839 * I) * q2 * q2) / (-0.06760000000000001 + q2)
                         + (20.10899364393711 + 17.72572790965085 * I) * std::log(mu)
                         - ((1.3498965971608112 + 0.22641335084230746 * I) * std::log(mu)) / (-0.06760000000000001 + q2)
                         + (3.9034896415281084 - 6.532993250678167 * I) * q2 * std::log(mu)
                         - ((14.797389926616242 + 0.3000399789915488 * I) * q2 * std::log(mu)) / (-0.06760000000000001 + q2)
                         + (0.3629832585569904 - 0.014240457252767418 * I) * q2 * q2 * std::log(mu)
                         - ((1.6426778147318113 + 6.4515613351795675 * I) * q2 * q2 * std::log(mu)) / (-0.06760000000000001 + q2)
                         + 0.45691388417234213 * std::log(mu) * std::log(mu) - (0.0024952022807061647 * std::log(mu) * std::log(mu)) / (-0.06760000000000001 + q2)
                         - 0.006791331988156049 * q2 * std::log(mu) * std::log(mu) + (0.03645567687923401 * q2 * std::log(mu) * std::log(mu)) / (-0.06760000000000001 + q2)
                         + (0.006789662824910638 * q2 * q2 * std::log(mu) * std::log(mu)) / (-0.06760000000000001 + q2)
                         + (38.740351545262385 + 51.323306973054116 * I) * std::log(-0.06760000000000001 + q2_c)
                         + ((0.021426489176010016 - 0.003993656630197696 * I) * std::log(-0.06760000000000001 + q2_c)) / (-0.06760000000000001 + q2)
                         + (54.51846051805137 - 197.9175134739021 * I) * q2 * std::log(-0.06760000000000001 + q2_c)
                         - (0.403777417080432 - 0.3474916863453839 * I) * q2 * q2 * std::log(-0.06760000000000001 + q2_c)
                         + (0.35063623218243956 + 11.597796340319157 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c)
                         + ((0.0039428769367046606 + 0.0017625604833091313 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c)) / (-0.06760000000000001 + q2)
                         + (25.91431005671744 - 20.712696820244503 * I) * q2 * std::log(mu) * std::log(-0.06760000000000001 + q2_c)
                         - (0.06396371059116047 - 0.02137270385715039 * I) * q2 * q2 * std::log(mu) * std::log(-0.06760000000000001 + q2_c)
                         - (0.0276630201207335 - 3.3712624295957982 * I) * std::log(-0.06760000000000001 + q2_c) * std::log(-0.06760000000000001 + q2_c)
                         - (0.3083662911036639 - 0.6549968350987246 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(-0.06760000000000001 + q2_c)
                         + (3.154744014717021 - 3.4206781291862716 * I) * std::log(q2) - ((3.384174432628124 - 0.770604646916672 * I) * std::log(q2)) / (-0.06760000000000001 + q2)
                         - (39.2120435038448 - 236.02568382481866 * I) * q2 * std::log(q2)
                         - ((47.79266039590835 - 10.472063873087809 * I) * q2 * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (0.9260865595069737 - 1.4644714975732458 * I) * std::log(mu) * std::log(q2)
                         - ((0.5094928428490364 + 0.20675518830339665 * I) * std::log(mu) * std::log(q2)) / (-0.06760000000000001 + q2)
                         - (26.39589328095337 - 26.252239392871395 * I) * q2 * std::log(mu) * std::log(q2)
                         - ((7.126693542734449 + 3.1717566696042114 * I) * q2 * std::log(mu) * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (8.11094659466384 + 7.441435333789422 * I) * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         - (4.920071189060938 + 4.424216267625104 * I) * q2 * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         + (0.5475454453374203 + 1.480981591645106 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         - (0.3317956261695305 + 0.8087472665751683 * I) * q2 * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         - (4.133715720547897 - 1.1692479729029264 * I) * std::log(q2) * std::log(q2)
                         - ((0.6247853354673722 + 0.4600024407081391 * I) * std::log(q2) * std::log(q2)) / (-0.06760000000000001 + q2)
                         - (0.5674513877714784 + 0.3538103275221192 * I) * std::log(mu) * std::log(q2) * std::log(q2)
                         - ((0.038509353091227674 + 0.11726786644229804 * I) * std::log(mu) * std::log(q2) * std::log(q2)) / (-0.06760000000000001 + q2)
                         - (0.8567928895868717 - 2.1633608294660576 * I) * std::log(-0.06760000000000001 + q2_c) * std::log(q2) * std::log(q2)
                         - (0.27582712466640935 - 0.19545708358270839 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(q2) * std::log(q2));
                complex<double> F_91_d =
                        ((-0.7028772801704175 - 1.1231469897640516 * I) + (0.015398873010909985 - 0.16381373572467825 * I) / (q2 * q2)
                         + (0.5417530556070507 + 1.5862030556023665 * I) / q2 + (0.6473419983308646 - 0.3191829251462698 * I) * q2
                         - (0.025313311807756117 - 0.05786601882115685 * I) * q2 * q2 - (0.9932540647119045 + 0.7235955663584045 * I) * std::log(mu)
                         + ((0.0048443954694649995 - 0.0034337873967769505 * I) * std::log(mu)) / q2 + (0.0023518768894868042 - 0.002329531067035292 * I) * q2 * std::log(mu)
                         + 0.0003815791570362812 * q2 * q2 * std::log(mu) - 0.0823043594593947 * std::log(mu) * std::log(mu)
                         + (1.3280134271551518 - 0.23019162732422474 * I) * std::log(-q2_c) + ((0.0003450321424401025 - 0.02253983434291395 * I) * std::log(-q2_c)) / (q2 * q2)
                         - ((0.7401833027636743 - 0.3036044910255162 * I) * std::log(-q2_c)) / q2 - (0.304310135036402 + 0.03460707127212294 * I) * q2 * std::log(-q2_c)
                         - 0.0025524352607723534 * std::log(mu) * std::log(-q2_c) + ((0.001296643556525919 - 0.0002052768267481986 * I) * std::log(mu) * std::log(-q2_c)) / q2
                         + (0.0007413691267320502 + 0.0002384066015513561 * I) * q2 * std::log(mu) * std::log(-q2_c)
                         - (0.04423486138401016 - 0.4027183167128223 * I) * std::log(-q2_c) * std::log(-q2_c)
                         + ((0.0005008722138681334 - 0.001687323191702575 * I) * std::log(-q2_c) * std::log(-q2_c)) / (q2 * q2)
                         - ((0.12566547325129496 + 0.021286535165585697 * I) * std::log(-q2_c) * std::log(-q2_c)) / q2
                         + (0.00024326371557747804 - 0.0007728426001637217 * I) * std::log(mu) * std::log(-q2_c) * std::log(-q2_c)
                         + (0.00018070642613636175 * std::log(mu) * std::log(-q2_c) * std::log(-q2_c)) / q2 + (2.366787746962223 + 1.2872879171622464 * I) * std::log(q2)
                         + ((0.01271329129843619 - 0.0220132096090156 * I) * std::log(q2)) / (q2 * q2) + ((0.8778602416246676 + 0.3791771019545544 * I) * std::log(q2)) / q2
                         + (0.17788658833295912 + 0.20404503514118752 * I) * q2 * std::log(q2) + (0.2327566821583774 - 0.0015968209157409925 * I) * std::log(mu) * std::log(q2)
                         + ((0.0008896093097014679 - 0.0009302539448979172 * I) * std::log(mu) * std::log(q2)) / q2
                         + (0.0007415551944040773 - 0.0002383521501199373 * I) * q2 * std::log(mu) * std::log(q2)
                         - (0.4601015018931022 - 0.25891406567857644 * I) * std::log(q2) * std::log(q2)
                         + ((0.0013444199089192024 - 0.001457769697786996 * I) * std::log(q2) * std::log(q2)) / (q2 * q2)
                         + ((0.09750403588592693 + 0.12041365494097124 * I) * std::log(q2) * std::log(q2)) / q2
                         + (0.0002431873087111308 + 0.0007730864197508181 * I) * std::log(mu) * std::log(q2) * std::log(q2)
                         + (0.0001807832593110863 * std::log(mu) * std::log(q2) * std::log(q2)) / q2);
                complex<double> F_91_s =
                        ((124.71946699856203 - 76.13392744656795 * I) - (3.868037572543674 - 7.67641251434955 * I) / (-0.06760000000000001 + q2)
                         - (35.00210181936058 + 38.23030469217163 * I) * q2 - ((32.20087409952672 - 89.95776448806568 * I) * q2) / (-0.06760000000000001 + q2)
                         + (0.5424635948512242 - 2.5861162747652595 * I) * q2 * q2 - ((45.795762096754146 + 4.08831101714455 * I) * q2 * q2) / (-0.06760000000000001 + q2)
                         - (3.7627207660727935 + 2.9542917020578434 * I) * std::log(mu)
                         + ((0.22722799745097655 + 0.037735663257519425 * I) * std::log(mu)) / (-0.06760000000000001 + q2)
                         - (0.6444651603174065 - 1.0888329081180688 * I) * q2 * std::log(mu)
                         + ((2.4334156847535158 + 0.05000748402911039 * I) * q2 * std::log(mu)) / (-0.06760000000000001 + q2)
                         - (0.060497140724599066 - 0.002373922268966332 * I) * q2 * q2 * std::log(mu)
                         + ((0.2676709112794735 + 1.075261202999089 * I) * q2 * q2 * std::log(mu)) / (-0.06760000000000001 + q2) - 0.07615234659821937 * std::log(mu) * std::log(mu)
                         + (0.00041588924556527137 * std::log(mu) * std::log(mu)) / (-0.06760000000000001 + q2) + 0.0011316817864188721 * q2 * std::log(mu) * std::log(mu)
                         - (0.006075657487386759 * q2 * std::log(mu) * std::log(mu)) / (-0.06760000000000001 + q2)
                         - (0.0011317053240523125 * q2 * q2 * std::log(mu) * std::log(mu)) / (-0.06760000000000001 + q2)
                         + (71.06988195712674 + 13.75714233373244 * I) * std::log(-0.06760000000000001 + q2_c)
                         + ((0.015413605587762846 - 0.019137455203434763 * I) * std::log(-0.06760000000000001 + q2_c)) / (-0.06760000000000001 + q2)
                         - (80.59545534225576 + 197.42375610798337 * I) * q2 * std::log(-0.06760000000000001 + q2_c)
                         + (0.04694994104091214 + 0.4890404358717765 * I) * q2 * q2 * std::log(-0.06760000000000001 + q2_c)
                         - (0.05843764769312331 + 1.9329652203815222 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c)
                         - ((0.0006571453385195534 + 0.0002937607377870104 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c)) / (-0.06760000000000001 + q2)
                         - (4.319048042608683 - 3.4521088215314077 * I) * q2 * std::log(mu) * std::log(-0.06760000000000001 + q2_c)
                         + (0.010660570753941872 - 0.0035622089451749377 * I) * q2 * q2 * std::log(mu) * std::log(-0.06760000000000001 + q2_c)
                         + (2.510596460700943 + 2.186545794762968 * I) * std::log(-0.06760000000000001 + q2_c) * std::log(-0.06760000000000001 + q2_c)
                         + (0.05139446163669032 - 0.10916603659890371 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(-0.06760000000000001 + q2_c)
                         - (1.07157361802855 + 6.83677159238133 * I) * std::log(q2) - ((2.2706380959072257 - 2.6716830993183387 * I) * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (116.47575004437326 + 212.18346280516116 * I) * q2 * std::log(q2)
                         - ((32.55143854938306 - 37.0699972473221 * I) * q2 * std::log(q2)) / (-0.06760000000000001 + q2)
                         - (0.15435046470187597 - 0.24407667052792734 * I) * std::log(mu) * std::log(q2)
                         + ((0.08491524278346001 + 0.03445918118865007 * I) * std::log(mu) * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (4.3993080338145205 - 4.375366744974276 * I) * q2 * std::log(mu) * std::log(q2)
                         + ((1.1877785506705296 + 0.5286255804163315 * I) * q2 * std::log(mu) * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (11.228714857886496 + 0.8994537069897728 * I) * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         - (6.072673749743096 - 0.6342293466615058 * I) * q2 * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         - (0.09125745419044744 + 0.24683046510344905 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         + (0.05529996628010756 + 0.13479089925273752 * I) * q2 * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(q2)
                         - (2.8255164988781423 - 3.2175128187359454 * I) * std::log(q2) * std::log(q2)
                         - ((0.7859261991845948 - 0.005340908557108383 * I) * std::log(q2) * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (0.09457452464607988 + 0.05896783490030841 * I) * std::log(mu) * std::log(q2) * std::log(q2)
                         + ((0.006418191263290735 + 0.01954460125379662 * I) * std::log(mu) * std::log(q2) * std::log(q2)) / (-0.06760000000000001 + q2)
                         + (0.6269131396060723 + 2.2068775272537255 * I) * std::log(-0.06760000000000001 + q2_c) * std::log(q2) * std::log(q2)
                         + (0.0459711223250179 - 0.032576212833344345 * I) * std::log(mu) * std::log(-0.06760000000000001 + q2_c) * std::log(q2) * std::log(q2));

                complex<double> C9eff_s = ret + c._alpha_s / (4 * M_PI) * Lms * (-8. / 27. * c.c1_i(1) - 2. / 9. * c.c2_i(1))
                                          + (c._alpha_s * c._alpha_s / (4 * M_PI * 4 * M_PI)) * (F_91_s * c.c1_i(0) + F_92_s * c.c2_i(0) + F_89 * C8_eff);
                complex<double> C9eff_d = ret + c._alpha_s / (4 * M_PI) * Lzero * (-8. / 27. * c.c1_i(1) - 2. / 9. * c.c2_i(1))
                                          + (c._alpha_s * c._alpha_s / (4 * M_PI * 4 * M_PI)) * (F_91_d * c.c1_i(0) + F_92_d * c.c2_i(0) + F_89 * C8_eff);

                const double V_cd       = abs_Vcd;
                const double V_cd_theta = arg_Vcd;
                const double V_ud       = abs_Vud;

                const double V_cs       = abs_Vcs;
                const double V_cs_theta = arg_Vcs;
                const double V_us       = abs_Vus;

                if (opt_cp_conjugate.value())
                {
                    return 4 * M_PI / c._alpha_s * (std::polar(V_cd, V_cd_theta) * V_ud * C9eff_d + std::polar(V_cs, V_cs_theta) * V_us * C9eff_s);
                }
                return 4 * M_PI / c._alpha_s * (std::polar(V_cd, -V_cd_theta) * V_ud * C9eff_d + std::polar(V_cs, -V_cs_theta) * V_us * C9eff_s);
            }

            lambda_c_to_proton_l_l::Amplitudes
            amplitudes(const double & s)
            {
                using std::conj;
                using std::imag;
                using std::norm;
                using std::real;
                using std::sqrt;

                lambda_c_to_proton_l_l::Amplitudes result;

                const auto wc    = model->wilson_coefficients_uc(opt_l.value(), opt_cp_conjugate.value());
                const auto wc_SM = model_SM->wilson_coefficients_uc(opt_l.value(), opt_cp_conjugate.value());

                const complex<double> C9_eff = C9eff(s, wc_SM);

                const complex<double> c7minus  = wc.c7() - wc.c7prime();
                const complex<double> c7plus   = wc.c7() + wc.c7prime();
                const complex<double> c9minus  = C9_eff + wc.c9() - wc.c9prime();
                const complex<double> c9plus   = C9_eff + wc.c9() + wc.c9prime();
                const complex<double> c10minus = wc.c10() - wc.c10prime();
                const complex<double> c10plus  = wc.c10() + wc.c10prime();

                const double res_a_rho   = a_rho;
                const double res_a_omega = a_omega;
                const double res_a_phi   = a_phi;

                const double res_delta_rho   = delta_rho;
                const double res_delta_omega = delta_omega_m_rho;
                const double res_delta_phi   = delta_phi_m_rho;

                const complex<double> c9R =
                        (res_a_rho / (complex<double>(s - m_rho * m_rho, m_rho * gamma_rho))
                         + res_a_omega * complex<double>(std::cos(res_delta_omega), std::sin(res_delta_omega)) / (complex<double>(s - m_omega * m_omega, m_omega * gamma_omega))
                         + res_a_phi * complex<double>(std::cos(res_delta_phi), std::sin(res_delta_phi)) / (complex<double>(s - m_phi * m_phi, m_phi * gamma_phi)));
                const complex<double> cPhase(std::cos(res_delta_rho), std::sin(res_delta_rho));

                // baryonic form factors (10)
                const double fff0         = form_factors->f_time_v(s);  // f0
                const double fffplus      = form_factors->f_long_v(s);  // fplus
                const double fffperp      = form_factors->f_perp_v(s);  // fperp
                const double ffg0         = form_factors->f_time_a(s);  // g0
                const double ffgplus      = form_factors->f_long_a(s);  // gplus
                const double ffgperp      = form_factors->f_perp_a(s);  // gperp
                const double ffhplus      = form_factors->f_long_t(s);  // hplus
                const double ffhtildeplus = form_factors->f_long_t5(s); // htildeplus
                const double ffhperp      = form_factors->f_perp_t(s);  // hperp
                const double ffhtildeperp = form_factors->f_perp_t5(s); // htildeperp
                // running quark masses
                const double mcatmu       = model->m_c_msbar(mu);

                // kinematics
                const double beta_l_squared   = (1.0 - 4.0 * m_l * m_l / s);
                const double beta_l           = std::sqrt(beta_l_squared);
                const double four_m_l_squared = 4.0 * m_l * m_l;
                const double mplus            = m_Lambda_c + m_proton;
                const double mminus           = m_Lambda_c - m_proton;
                const double mplus_squared    = power_of<2>(mplus);
                const double mminus_squared   = power_of<2>(mminus);
                const double splus            = mplus_squared - s;
                const double sminus           = mminus_squared - s;
                const double sqrtsminussplus  = std::sqrt(sminus) * std::sqrt(splus);

                // normalization
                const complex<double> N = normalization(s);

                // cf. [GHM:2021A], eq. (11), p. 6, excluding contributions from scalar and tensor operators.
                result.aU11 = +4.0
                              * ((norm(c7plus * 2.0 * mcatmu / s * mplus * ffhperp + c9plus * fffperp)
                                  + 2.0 * real(conj(c7plus * 2.0 * mcatmu / s * mplus * ffhperp + c9plus * fffperp) * c9R * cPhase * fffperp) + norm(c9R * fffperp))
                                         * sminus
                                 + (norm(c7minus * 2.0 * mcatmu / s * mminus * ffhtildeperp + c9minus * ffgperp)
                                    + 2.0 * real(conj(c7minus * 2.0 * mcatmu / s * mminus * ffhtildeperp + c9minus * ffgperp) * c9R * cPhase * ffgperp) + norm(c9R * ffgperp))
                                           * splus);

                result.aL11 = +2.0 / s
                              * ((norm(c7plus * 2.0 * mcatmu * ffhplus + c9plus * mplus * fffplus)
                                  + 2.0 * real(conj(c7plus * 2.0 * mcatmu * ffhplus + c9plus * mplus * fffplus) * c9R * cPhase * mplus * fffplus) + norm(c9R * mplus * fffplus))
                                         * sminus
                                 + (norm(c7minus * 2.0 * mcatmu * ffhtildeplus + c9minus * mminus * ffgplus)
                                    + 2.0 * real(conj(c7minus * 2.0 * mcatmu * ffhtildeplus + c9minus * mminus * ffgplus) * c9R * cPhase * mminus * ffgplus)
                                    + norm(c9R * mminus * ffgplus))
                                           * splus);
                result.aU22 = +4.0 * (norm(c10plus) * fffperp * fffperp * sminus + norm(c10minus) * ffgperp * ffgperp * splus);
                result.aL22 = +2.0 / s * (norm(c10plus) * mplus_squared * fffplus * fffplus * sminus + norm(c10minus) * mminus_squared * ffgplus * ffgplus * splus);
                result.aS22 = +2.0 / s * (norm(c10plus) * mminus_squared * fff0 * fff0 * splus + norm(c10minus) * mplus_squared * ffg0 * ffg0 * sminus);
                result.aP12 =
                        -8
                        * (real(c7minus * conj(c10plus)) * mcatmu / s * mminus * fffperp * ffhtildeperp + real(c7plus * conj(c10minus)) * mcatmu / s * mplus * ffgperp * ffhperp
                           + real((wc.c9() + c9R * cPhase) * conj(wc.c10()) - wc.c9prime() * conj(wc.c10prime())) * ffgperp * fffperp)
                        * sqrtsminussplus;

                result.beta_l           = beta_l;
                result.beta_l_squared   = beta_l_squared;
                result.four_m_l_squared = four_m_l_squared;
                result.q2               = s;
                result.norm_squared     = norm(N);

                return result;
            }

            std::array<double, 3>
            _differential_angular_observables(const double & q2)
            {
                return lambda_c_to_proton_l_l::AngularObservables(this->amplitudes(q2))._k;
            }

            std::array<double, 3>
            _integrated_angular_observables(const double & q2_min, const double & q2_max)
            {
                std::function<std::array<double, 3>(const double &)> integrand = [this](const double & q2) -> std::array<double, 3>
                { return this->_differential_angular_observables(q2); };

                return integrate<1, 3>(integrand, q2_min, q2_max, cubature::Config().epsrel(1e-5));
                ;
            }

            inline lambda_c_to_proton_l_l::AngularObservables
            differential_angular_observables(const double & q2)
            {
                return lambda_c_to_proton_l_l::AngularObservables{ _differential_angular_observables(q2) };
            }

            inline lambda_c_to_proton_l_l::AngularObservables
            integrated_angular_observables(const double & q2_min, const double & q2_max)
            {
                return lambda_c_to_proton_l_l::AngularObservables{ _integrated_angular_observables(q2_min, q2_max) };
            }
    };

    const std::vector<OptionSpecification> Implementation<LambdaCToProtonLeptonLepton>::options{
        Model::option_specification(),
        { "cp-conjugate"_ok,   { "true"s, "false"s }, "false"s },
        {            "l"_ok, { "e"s, "mu"s, "tau"s },    "mu"s },
    };

    LambdaCToProtonLeptonLepton::LambdaCToProtonLeptonLepton(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<LambdaCToProtonLeptonLepton>(new Implementation<LambdaCToProtonLeptonLepton>(p, o, *this))
    {
    }

    LambdaCToProtonLeptonLepton::~LambdaCToProtonLeptonLepton() {}

    /* for double-differential signal PDF */
    double
    LambdaCToProtonLeptonLepton::double_differential_decay_width(const double & q2, const double & c_lep) const
    {
        return _imp->differential_angular_observables(q2).d2gamma(c_lep);
    }

    double
    LambdaCToProtonLeptonLepton::integrated_decay_width(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).decay_width();
    }

    /* q^2-differential observables */

    double
    LambdaCToProtonLeptonLepton::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToProtonLeptonLepton::differential_a_fb_leptonic(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_leptonic();
    }

    double
    LambdaCToProtonLeptonLepton::differential_f_l(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).f_l();
    }

    /* q^2-integrated observables */

    double
    LambdaCToProtonLeptonLepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic_num(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic_num();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_f_l(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_l();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_f_l_num(const double & s_min, const double & s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_l_num();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1ss(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1ss() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1cc(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1cc() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1c(const double & s_min, const double & s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1c() / o.decay_width();
    }

    const std::string LambdaCToProtonLeptonLepton::description = "\
      The decay Lambda_c -> proton lbar l, where lbar=e^+,mu^+,tau^+ is a charged antilepton.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_q2 = "\
      The invariant mass of the lbar-l pair in GeV^2.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_c_theta_l = "\
      The cosine of the helicity angle between the direction of flight of the charged antilepton and of the Lambda_c in the lbar-l rest frame.";

    const std::set<ReferenceName> LambdaCToProtonLeptonLepton::references{};

    std::vector<OptionSpecification>::const_iterator
    LambdaCToProtonLeptonLepton::begin_options()
    {
        return Implementation<LambdaCToProtonLeptonLepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    LambdaCToProtonLeptonLepton::end_options()
    {
        return Implementation<LambdaCToProtonLeptonLepton>::options.cend();
    }
} // namespace eos

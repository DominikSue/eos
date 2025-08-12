/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Ahmet Kokulu
 * Copyright (c) 2019,2021 Danny van Dyk
 * Copyright (c) 2023 Méril Reboud
 * Copyright (c) 2025 Dominik Suelmann
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

#include <eos/c-decays/lambdac-to-proton-l-l.hh>
#include <eos/form-factors/baryonic.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
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
    using std::norm;

    namespace lambdac_to_proton_l_l
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
            std::array<double, 10> _k;

            AngularObservables(const Amplitudes &a)
            {
                using std::conj;
                using std::imag;
                using std::norm;
                using std::real;
                using std::sqrt;

                double beta_l = a.beta_l;
                double beta_l_squared = a.beta_l_squared;
                double four_m_l_squared = a.four_m_l_squared;
                double q2 = a.q2;
                double aU11U22 = a.aU11 + a.aU22;
                double aULS = a.aU11 + a.aL11 + a.aS22;

                // cf. [BKTvD2019], eqs. (2.5)-(2.14), pp. 4-5.

                // K_{1ss}
                _k[0] = a.norm_squared * (q2 * beta_l_squared * (aU11U22 / 2.0 + a.aL11 + a.aL22) + four_m_l_squared * aULS);

                // K_{1cc}
                _k[1] = a.norm_squared * (q2 * beta_l_squared * aU11U22 + four_m_l_squared * aULS);

                // K_{1c}
                _k[2] = a.norm_squared * (-2.0 * q2 * beta_l * a.aP12);

                // K_{11}
                _k[3] = 0.0;

                // K_{12}
                _k[4] = 0.0;

                // K_{13}
                _k[5] = 0.0;

                // K_{21}
                _k[6] = 0.0;

                // K_{23}
                _k[7] = 0.0;

                // K_{22}
                _k[8] = 0.0;

                // K_{24}
                _k[9] = 0.0;
            }

            AngularObservables(const std::array<double, 10> &k) : _k(k)
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
            k2ss() const
            {
                return _k[3];
            }

            inline double
            k2cc() const
            {
                return _k[4];
            }

            inline double
            k2c() const
            {
                return _k[5];
            }

            inline double
            k3sc() const
            {
                return _k[6];
            }

            inline double
            k3s() const
            {
                return _k[7];
            }

            inline double
            k4sc() const
            {
                return _k[8];
            }

            inline double
            k4s() const
            {
                return _k[9];
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
            a_fb_hadronic() const
            {
                return 1.0 / 2.0 * (2.0 * k2ss() + k2cc()) / decay_width();
            }

            inline double
            a_fb_combined() const
            {
                return 3.0 / 4.0 * k2c() / decay_width();
            }

            inline double
            f_zero() const
            {
                return (2.0 * k1ss() - k1cc()) / decay_width();
            }

            inline double
            f_L_num() const
            {
                return (2.0 * k1ss() - k1cc());
            }

            inline double
            d4gamma(const double &c_lep, const double &c_lam, const double &phi) const
            {
                const double c2_lep = c_lep * c_lep;
                const double s2_lep = 1.0 - c2_lep;
                const double s_lep = std::sqrt(s2_lep);
                const double s_lam = std::sqrt(1.0 - c_lam * c_lam);
                const double c_phi = std::cos(phi), s_phi = std::sin(phi);

                // cf. [BKTvD2019], p. 2, eqs. (2.3) and (2.4)
                return 3.0 / (8.0 * M_PI) * (k1ss() * s2_lep + k1cc() * c2_lep + k1c() * c_lep + (k2ss() * s2_lep + k2cc() * c2_lep + k2c() * c_lep) * c_lam + (k3sc() * s_lep * c_lep + k3s() * s_lep) * s_lam * s_phi + (k4sc() * s_lep * c_lep + k4s() * s_lep) * s_lam * c_phi);
            }
        };
    } // namespace lambdac_to_proton_l_l

    template <>
    struct Implementation<LambdaCToProtonLeptonLepton>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        UsedParameter hbar;
        UsedParameter tau_Lambda_c;

        UsedParameter g_fermi;
        UsedParameter alpha_e;

        LeptonFlavorOption opt_l;
        UsedParameter m_l;

        UsedParameter m_Lambda_c;
        UsedParameter m_Proton;
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

        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> form_factors;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters &p, const Options &o, ParameterUser &u) : model(Model::make(o.get("model"_ok, "SM"), p, o)),
                                                                                  parameters(p),
                                                                                  hbar(p["QM::hbar"], u),
                                                                                  tau_Lambda_c(p["life_time::Lambda_c"], u),
                                                                                  g_fermi(p["WET::G_Fermi"], u),
                                                                                  alpha_e(p["QED::alpha_e(m_c)"], u),
                                                                                  opt_l(o, options, "l"_ok),
                                                                                  m_l(p["mass::" + opt_l.str()], u),
                                                                                  m_Lambda_c(p["mass::Lambda_c"], u),
                                                                                  m_Proton(p["mass::Proton"], u),
                                                                                  m_rho(p["mass::rho^0"], u),
                                                                                  m_omega(p["mass::omega"], u),
                                                                                  m_phi(p["mass::phi"], u),
                                                                                  gamma_rho(p["gamma::rho^0"], u),
                                                                                  gamma_omega(p["gamma::omega"], u),
                                                                                  gamma_phi(p["gamma::phi"], u),
                                                                                  a_rho(p["Lambda_c->Proton::res_a_rho"], u),
                                                                                  a_omega(p["Lambda_c->Proton::res_a_omega"], u),
                                                                                  a_phi(p["Lambda_c->Proton::res_a_phi"], u),
                                                                                  delta_rho(p["Lambda_c->Proton::res_delta_rho"], u),
                                                                                  delta_omega_m_rho(p["Lambda_c->Proton::res_delta_omega_m_rho"], u),
                                                                                  delta_phi_m_rho(p["Lambda_c->Proton::res_delta_phi_m_rho"], u),
                                                                                  opt_cp_conjugate(o, options, "cp-conjugate"_ok),
                                                                                  mu(p["uc" + opt_l.str() + opt_l.str() + "::mu"], u),
                                                                                  form_factors(FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_c->Proton::" + o.get("form-factors"_ok, "DM2016"), p, o))
        {
            u.uses(*form_factors);
            u.uses(*model);
        }

        const complex<double>
        normalization(const double &q2) const
        {
            double lam = lambda(m_Lambda_c * m_Lambda_c, m_Proton * m_Proton, q2);

            if ((lam <= 0) || (q2 <= 4.0 * m_l * m_l))
            {
                return 0.0;
            }

            // std::printf("alpha_e: %f\n", alpha_e());
            // std::printf("g_fermi: %e\n", g_fermi());

            // alpha_e = 0.007501720229738859
            return g_fermi() * 0.007501720229738859 * std::sqrt(std::sqrt(1.0 - 4.0 * m_l * m_l / q2)) * std::sqrt(std::sqrt(lam) / 3.0 / power_of<11>(2) / power_of<5>(M_PI) / power_of<3>(m_Lambda_c()));
        }

        lambdac_to_proton_l_l::Amplitudes
        amplitudes(const double &s)
        {
            using std::conj;
            using std::imag;
            using std::norm;
            using std::real;
            using std::sqrt;

            lambdac_to_proton_l_l::Amplitudes result;

            const auto wc = model->wet_ucll(opt_l.value(), opt_cp_conjugate.value()); // need to redo
            const complex<double> c7 = wc.c7();
            const complex<double> c7p = wc.c7p();
            const complex<double> c9 = wc.c9();
            const complex<double> c9p = wc.c9p();
            const complex<double> c10 = wc.c10();
            const complex<double> c10p = wc.c10p();

            const complex<double> c7minus = c7-c7p;
            const complex<double> c7plus = c7+c7p;
            const complex<double> c9minus = c9-c9p;
            const complex<double> c9plus = c9+c9p;
            const complex<double> c10minus = c10-c10p;
            const complex<double> c10plus = c10+c10p;

            const complex<double> rho_fac = 1.0 / (s - m_rho * m_rho + std::polar(m_rho * gamma_rho, M_PI / 2.0));
            const complex<double> omega_fac = 1.0 / (s - m_omega * m_omega + std::polar(m_omega * gamma_omega, M_PI / 2.0));
            const complex<double> phi_fac = 1.0 / (s - m_phi * m_phi + std::polar(m_phi * gamma_phi, M_PI / 2.0));

            const double res_a_rho = a_rho;
            const double res_a_omega = a_omega;
            const double res_a_phi = a_phi;

            const double res_delta_rho = delta_rho;
            const double res_delta_omega = delta_omega_m_rho;
            const double res_delta_phi = delta_phi_m_rho;

            const complex<double> c9R = (res_a_rho * rho_fac + std::polar(res_a_omega, res_delta_omega) * omega_fac + std::polar(res_a_phi, res_delta_phi) * phi_fac);

            // std::printf("a_rho = %f, a_omega = %f, a_phi = %f\n", res_a_rho, res_a_omega, res_a_phi);
            // std::printf("delta_rho = %f, delta_omega = %f, delta_phi = %f\n", res_delta_rho, res_delta_omega, res_delta_phi);
            // std::printf("m_rho = %f, m_omega = %f, m_phi = %f\n", m_rho(), m_omega(), m_phi());
            // std::printf("gamma_rho = %f, gamma_omega = %f, gamma_phi = %f\n", gamma_rho(), gamma_omega(), gamma_phi());
            // std::printf("C9R(q2:%f) = %f + %fi\n", s, real(c9R), imag(c9R));

            // std::printf("rho_fac(q2:%f) = %f + %fi\n", s, real(rho_fac), imag(rho_fac));
            // std::printf("omega_fac(q2:%f) = %f + %fi\n", s, real(omega_fac), imag(omega_fac));
            // std::printf("phi_fac(q2:%f) = %f + %fi\n", s, real(phi_fac), imag(phi_fac));

            // std::printf("c7 = %f + %fi\n", real(c7), imag(c7));
            // std::printf("c7p = %f + %fi\n", real(c7p), imag(c7p));
            // std::printf("c9 = %f + %fi\n", real(c9), imag(c9));
            // std::printf("c9p = %f + %fi\n", real(c9p), imag(c9p));
            // std::printf("c10 = %f + %fi\n", real(c10), imag(c10));
            // std::printf("c10p = %f + %fi\n", real(c10p), imag(c10p));

            const complex<double> cPhase = std::polar(1.0, res_delta_rho);

            // std::printf("cPhase:%f + %fi\n", real(cPhase), imag(cPhase));

            // baryonic form factors (10)
            const double fff0 = form_factors->f_time_v(s);          // f0
            const double fffplus = form_factors->f_long_v(s);       // fplus
            const double fffperp = form_factors->f_perp_v(s);       // fperp
            const double ffg0 = form_factors->f_time_a(s);          // g0
            const double ffgplus = form_factors->f_long_a(s);       // gplus
            const double ffgperp = form_factors->f_perp_a(s);       // gperp
            const double ffhplus = form_factors->f_long_t(s);       // hplus
            const double ffhtildeplus = form_factors->f_long_t5(s); // htildeplus
            const double ffhperp = form_factors->f_perp_t(s);       // hperp
            const double ffhtildeperp = form_factors->f_perp_t5(s); // htildeperp
            // running quark masses
            const double mcatmu = model->m_c_msbar(mu);
            // const double msatmu = model->m_s_msbar(mu);

            // kinematics
            const double beta_l_squared = (1.0 - 4.0 * m_l * m_l / s);
            const double beta_l = std::sqrt(beta_l_squared);
            const double four_m_l_squared = 4.0 * m_l * m_l;
            // const double sqrts      = std::sqrt(s);
            const double mplus = m_Lambda_c + m_Proton;
            const double mminus = m_Lambda_c - m_Proton;
            const double mplus_squared = power_of<2>(mplus);
            const double mminus_squared = power_of<2>(mminus);
            const double splus = mplus_squared - s;
            const double sminus = mminus_squared - s;
            const double sqrtsminussplus = std::sqrt(sminus) * std::sqrt(splus);

            // normalization
            const complex<double> N = normalization(s);

            // std::printf("c7plus = %f + %f\n", real(c7plus), imag(c7plus));
            // std::printf("c7minus = %f + %f\n", real(c7minus), imag(c7minus));
            // std::printf("c9plus = %f + %f\n", real(c9plus), imag(c9plus));
            // std::printf("c9minus = %f + %f\n", real(c9minus), imag(c9minus));

            // b->c case transversity amplitudes
            // cf. [BKTvD2019], eqs. (2.18)-(2.23), p. 6, including contributions from the vector and scalar operators.
            result.aU11 = +4.0 * ((norm(c7plus * 2.0 * mcatmu / s * mplus * ffhperp + c9plus * fffperp) + 2.0 * real(conj(c7plus * 2.0 * mcatmu / s * mplus * ffhperp + c9plus * fffperp) * c9R * cPhase * fffperp) + norm(c9R * fffperp)) * sminus + (norm(c7minus * 2.0 * mcatmu / s * mminus * ffhtildeperp + c9minus * ffgperp) + 2.0 * real(conj(c7minus * 2.0 * mcatmu / s * mminus * ffhtildeperp + c9minus * ffgperp) * c9R * cPhase * ffgperp) + norm(c9R * ffgperp)) * splus);

            result.aL11 = +2.0 / s * ((norm(c7plus * 2.0 * mcatmu * ffhplus + c9plus * mplus * fffplus) + 2.0 * real(conj(c7plus * 2.0 * mcatmu * ffhplus + c9plus * mplus * fffplus) * c9R * cPhase * mplus * fffplus) + norm(c9R * mplus * fffplus)) * sminus + (norm(c7minus * 2.0 * mcatmu * ffhtildeplus + c9minus * mminus * ffgplus) + 2.0 * real(conj(c7minus * 2.0 * mcatmu * ffhtildeplus + c9minus * mminus * ffgplus) * c9R * cPhase * mminus * ffgplus) + norm(c9R * mminus * ffgplus)) * splus);
            result.aU22 = +4.0 * (norm(c10plus) * fffperp * fffperp * sminus + norm(c10minus) * ffgperp * ffgperp * splus);
            result.aL22 = +2.0 / s * (norm(c10plus) * mplus_squared * fffplus * fffplus * sminus + norm(c10minus) * mminus_squared * ffgplus * ffgplus * splus);
            result.aS22 = +2.0 / s * (norm(c10plus) * mminus_squared * fff0 * fff0 * splus + norm(c10minus) * mplus_squared * ffg0 * ffg0 * sminus);
            result.aP12 =
                -8 * (real(c7minus * conj(c10plus)) * mcatmu / s * mminus * fffperp * ffhtildeperp + real(c7plus * conj(c10minus)) * mcatmu / s * mplus * ffgperp * ffhperp + real((c9 + c9R * cPhase) * conj(c10) - c9p * conj(c10p)) * ffgperp * fffperp) * sqrtsminussplus;

            // std::printf("fperp(q2=%f) = %f\n", s, fffperp);
            // std::printf("fplus(q2=%f) = %f\n", s, fffplus);
            // std::printf("ffhperp(q2=%f) = %f\n", s, ffhperp);
            // std::printf("ffgperp(q2=%f) = %f\n", s, ffgperp);
            // std::printf("ffhplus(q2=%f) = %f\n", s, ffhplus);
            // std::printf("ffhtildeperp(q2=%f) = %f\n", s, ffhtildeperp);
            // std::printf("ffgplus(q2=%f) = %f\n", s, ffgplus);
            // std::printf("mplus = %f, mminus = %f \n", mplus, mminus);

            result.beta_l = beta_l;
            result.beta_l_squared = beta_l_squared;
            result.four_m_l_squared = four_m_l_squared;
            result.q2 = s;
            result.norm_squared = norm(N);

            // std::printf("U11(q2=%f) = %f\n", s, result.aU11);
            // std::printf("L11(q2=%f) = %f\n", s, result.aL11);
            // std::printf("U22(q2=%f) = %f\n", s, result.aU22);
            // std::printf("L22(q2=%f) = %f\n", s, result.aL22);
            // std::printf("S22(q2=%f) = %f\n", s, result.aS22);
            // std::printf("P12(q2=%f) = %f\n", s, result.aP12);

            // std::printf("U11(q2=%f) = %e\n", s, result.aU11 * result.norm_squared);
            // std::printf("L11(q2=%f) = %e\n", s, result.aL11 * result.norm_squared);
            // std::printf("U22(q2=%f) = %e\n", s, result.aU22 * result.norm_squared);
            // std::printf("L22(q2=%f) = %e\n", s, result.aL22 * result.norm_squared);
            // std::printf("S22(q2=%f) = %e\n", s, result.aS22 * result.norm_squared);
            // std::printf("P12(q2=%f) = %e\n", s, result.aP12 * result.norm_squared);

            // std::printf("q2=%f, betal:%e, betal_squared:%e, 4ml^2:%e, |N|^2:%e \n", s, beta_l, beta_l_squared, four_m_l_squared, norm(N));
            return result;
        }

        std::array<double, 10>
        _differential_angular_observables(const double &q2)
        {
            return lambdac_to_proton_l_l::AngularObservables(this->amplitudes(q2))._k;
        }

        std::array<double, 10>
        _integrated_angular_observables(const double &q2_min, const double &q2_max)
        {
            std::function<std::array<double, 10>(const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));

            return integrate1D(integrand, 64, q2_min, q2_max);
        }

        inline lambdac_to_proton_l_l::AngularObservables
        differential_angular_observables(const double &q2)
        {
            return lambdac_to_proton_l_l::AngularObservables{_differential_angular_observables(q2)};
        }

        inline lambdac_to_proton_l_l::AngularObservables
        integrated_angular_observables(const double &q2_min, const double &q2_max)
        {
            return lambdac_to_proton_l_l::AngularObservables{_integrated_angular_observables(q2_min, q2_max)};
        }
    };

    const std::vector<OptionSpecification> Implementation<LambdaCToProtonLeptonLepton>::options{
        { "cp-conjugate"_ok, { "true", "false" },                     "false" },
        {"l"_ok, {"e", "mu", "tau"}, "mu"}};

    LambdaCToProtonLeptonLepton::LambdaCToProtonLeptonLepton(const Parameters &p, const Options &o) : PrivateImplementationPattern<LambdaCToProtonLeptonLepton>(new Implementation<LambdaCToProtonLeptonLepton>(p, o, *this))
    {
    }

    LambdaCToProtonLeptonLepton::~LambdaCToProtonLeptonLepton() {}

    /* for four-differential signal PDF */
    double
    LambdaCToProtonLeptonLepton::four_differential_decay_width(const double &q2, const double &c_lep, const double &c_lam, const double &phi) const
    {
        return _imp->differential_angular_observables(q2).d4gamma(c_lep, c_lam, phi);
    }

    double
    LambdaCToProtonLeptonLepton::integrated_decay_width(const double &q2_min, const double &q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).decay_width();
    }

    /* q^2-differential observables */

    double
    LambdaCToProtonLeptonLepton::differential_branching_ratio(const double &q2) const
    {
        // std::printf("tau_Lambda_c: %e\n", _imp->tau_Lambda_c / _imp->hbar);
        return _imp->differential_angular_observables(q2).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToProtonLeptonLepton::differential_a_fb_leptonic(const double &q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_leptonic();
    }

    double
    LambdaCToProtonLeptonLepton::differential_a_fb_hadronic(const double &q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_hadronic();
    }

    double
    LambdaCToProtonLeptonLepton::differential_a_fb_combined(const double &q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_combined();
    }

    double
    LambdaCToProtonLeptonLepton::differential_fzero(const double &q2) const
    {
        return _imp->differential_angular_observables(q2).f_zero();
    }

    /* q^2-integrated observables */

    double
    LambdaCToProtonLeptonLepton::integrated_branching_ratio(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).decay_width() * _imp->tau_Lambda_c / _imp->hbar;
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_leptonic_num(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_leptonic_num();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_hadronic(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_hadronic();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_a_fb_combined(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).a_fb_combined();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_fzero(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_zero();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_fL_num(const double &s_min, const double &s_max) const
    {
        return _imp->integrated_angular_observables(s_min, s_max).f_L_num();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1ss(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1ss() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1cc(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1cc() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k1c(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k1c() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k2ss(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2ss() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k2cc(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2cc() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k2c(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k2c() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k3sc(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3sc() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k3s(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k3s() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k4sc(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4sc() / o.decay_width();
    }

    double
    LambdaCToProtonLeptonLepton::integrated_k4s(const double &s_min, const double &s_max) const
    {
        auto o = _imp->integrated_angular_observables(s_min, s_max);
        return o.k4s() / o.decay_width();
    }

    const std::string LambdaCToProtonLeptonLepton::description = "\
      The decay Lambda_c -> Proton lbar l, where lbar=e^+,mu^+,tau^+ is a charged antilepton.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_q2 = "\
      The invariant mass of the lbar-l pair in GeV^2.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_c_theta_l = "\
      The cosine of the helicity angle between the direction of flight of the charged antilepton and of the Lambda_c in the lbar-l rest frame.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_c_theta_L = "\
      The cosine of the helicity angle between the direction of flight of the Lambda and of the pion in the Lambda_c rest frame.";

    const std::string LambdaCToProtonLeptonLepton::kinematics_description_phi = "\
      The azimuthal angle between the two decay planes.";

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

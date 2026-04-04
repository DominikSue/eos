/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/rare-c-decays/lambda-c-to-proton-l-l.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <test/test.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class LambdaCToProtonLeptonLeptonTest : public TestCase
{
    public:
        LambdaCToProtonLeptonLeptonTest() :
            TestCase("lambdac_to_proton_l_l_test")
        {
        }

        virtual void
        run() const
        {
            {
                Parameters p = Parameters::Defaults();

                p["Lambda_c->neutron::a^(t,V)_0@BMRvD2022"]     = +0.83558637190456;
                p["Lambda_c->neutron::a^(t,V)_1@BMRvD2022"]     = -2.5698251471323;
                p["Lambda_c->neutron::a^(t,V)_2@BMRvD2022"]     = +9.8729125261567;
                p["Lambda_c->neutron::a^(0,V)_0@BMRvD2022"]     = +0.83253055920653;
                p["Lambda_c->neutron::a^(0,V)_1@BMRvD2022"]     = -2.3309542260103;
                p["Lambda_c->neutron::a^(0,V)_2@BMRvD2022"]     = +8.4088373584501;
                p["Lambda_c->neutron::a^(perp,V)_0@BMRvD2022"]  = +1.3629101111921;
                p["Lambda_c->neutron::a^(perp,V)_1@BMRvD2022"]  = -1.6996176062383;
                p["Lambda_c->neutron::a^(perp,V)_2@BMRvD2022"]  = +0.70894108564687;
                p["Lambda_c->neutron::a^(t,A)_0@BMRvD2022"]     = +0.72755668638733;
                p["Lambda_c->neutron::a^(t,A)_1@BMRvD2022"]     = -0.96674724145661;
                p["Lambda_c->neutron::a^(t,A)_2@BMRvD2022"]     = +0.82644748408508;
                p["Lambda_c->neutron::a^(0,A)_0@BMRvD2022"]     = +0.68653490341925;
                p["Lambda_c->neutron::a^(0,A)_1@BMRvD2022"]     = -0.90303820170714;
                p["Lambda_c->neutron::a^(0,A)_2@BMRvD2022"]     = +2.249274327598;
                p["Lambda_c->neutron::a^(perp,A)_0@BMRvD2022"]  = +0.68653490341925;
                p["Lambda_c->neutron::a^(perp,A)_1@BMRvD2022"]  = -0.68345507203321;
                p["Lambda_c->neutron::a^(perp,A)_2@BMRvD2022"]  = +0.69589476713259;
                p["Lambda_c->neutron::a^(0,T)_0@BMRvD2022"]     = +1.111320720066;
                p["Lambda_c->neutron::a^(0,T)_1@BMRvD2022"]     = -0.68870038548072;
                p["Lambda_c->neutron::a^(0,T)_2@BMRvD2022"]     = -2.8436650308711;
                p["Lambda_c->neutron::a^(perp,T)_0@BMRvD2022"]  = +0.63353942203683;
                p["Lambda_c->neutron::a^(perp,T)_1@BMRvD2022"]  = -1.0352768906776;
                p["Lambda_c->neutron::a^(perp,T)_2@BMRvD2022"]  = +1.4209332847548;
                p["Lambda_c->neutron::a^(0,T5)_0@BMRvD2022"]    = +0.62562392794376;
                p["Lambda_c->neutron::a^(0,T5)_1@BMRvD2022"]    = -1.1924360577673;
                p["Lambda_c->neutron::a^(0,T5)_2@BMRvD2022"]    = +3.7297035832191;
                p["Lambda_c->neutron::a^(perp,T5)_0@BMRvD2022"] = +0.62562392794376;
                p["Lambda_c->neutron::a^(perp,T5)_1@BMRvD2022"] = -1.3886937865559;
                p["Lambda_c->neutron::a^(perp,T5)_2@BMRvD2022"] = +4.2242261969563;

                p["Lambda_c->proton::res_a_rho"]             = +0.54;
                p["Lambda_c->proton::res_a_omega"]           = +0.074;
                p["Lambda_c->proton::res_a_phi"]             = +0.106;
                p["Lambda_c->proton::res_delta_rho"]         = +0.00;
                p["Lambda_c->proton::res_delta_omega_m_rho"] = +M_PI;
                p["Lambda_c->proton::res_delta_phi_m_rho"]   = +M_PI;
                p["gamma::rho^0"]                            = 0.1474;

                Options oo{
                    {        "model"_ok,       "WET" },
                    { "form-factors"_ok, "BMRvD2022" },
                    {            "l"_ok,        "mu" }
                };

                LambdaCToProtonLeptonLepton d(p, oo);

                const double eps = 3e-2;

                TEST_CHECK_RELATIVE_ERROR(d.differential_branching_ratio(0.75), 2.75291755e-7, eps);

                // the full phase-space region for muon
                TEST_CHECK_RELATIVE_ERROR(d.integrated_branching_ratio(0.959, 1.122), 3.02e-7, eps);

                p["ucmumu::Re{c10}"] = +2.0;
                p["ucmumu::Im{c10}"] = 1.4;

                Options oobar{
                    {        "model"_ok,       "WET" },
                    { "form-factors"_ok, "BMRvD2022" },
                    {            "l"_ok,        "mu" },
                    { "cp-conjugate"_ok,      "true" }
                };

                LambdaCToProtonLeptonLepton dNP(p, oo);
                LambdaCToProtonLeptonLepton dNPbar(p, oobar);

                // check NP contributions
                TEST_CHECK_RELATIVE_ERROR(dNP.differential_a_fb_leptonic(1.5), 0.04848, eps);
                TEST_CHECK_RELATIVE_ERROR(dNP.differential_branching_ratio(1.5), 1.7868688553844158e-07, eps);
                TEST_CHECK_RELATIVE_ERROR(dNP.integrated_decay_width(0.4, 0.9), 9.309968605287797e-19, eps);
                // Sigma
                TEST_CHECK_RELATIVE_ERROR(0.5
                                                  * (dNP.integrated_a_fb_leptonic_num(0.4, 0.9) / dNP.integrated_decay_width(0.4, 0.9)
                                                     + dNPbar.integrated_a_fb_leptonic_num(0.4, 0.9) / dNPbar.integrated_decay_width(0.4, 0.9)),
                                          0.052767573604217216,
                                          eps);
                // Delta
                TEST_CHECK_RELATIVE_ERROR(0.5
                                                  * (dNP.integrated_a_fb_leptonic_num(0.4, 0.9) / dNP.integrated_decay_width(0.4, 0.9)
                                                     - dNPbar.integrated_a_fb_leptonic_num(0.4, 0.9) / dNPbar.integrated_decay_width(0.4, 0.9)),
                                          -0.0786251567266188,
                                          eps);
            }
        }
} lambdac_to_proton_l_l_test;

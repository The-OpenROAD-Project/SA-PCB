///////////////////////////////////////////////////////////////////////////////
// Authors: Chester Holtz, Devon Merrill, James (Ting-Chou) Lin, Connie (Yen-Yi) Wu
//          (respective Ph.D. advisors: Chung-Kuan Cheng, Andrew B. Kahng, Steven Swanson).
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#define BOOST_NO_AUTO_PTR

#include <cstdio>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <map>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <numeric>
#include <cmath>
#include <unistd.h>

class Logger {
	public:
		vector<double> acceptance_ratio_history;
		vector<double> temperature_history;
		vector<double> move_window_coef_history;
		vector<double> cost_history;
		vector<double> hpwl_history;
		vector<double> overlap_history;
		vector<double> overlap_x_history;
                vector<double> overlap_x2_history;
	        vector<double> cdist_history;	
		vector<double> congestion_history;
		vector<double> lambda_history;
		vector<double> lamrate_history;

		string micro_name ="c31";
		vector<double> micro_dcost;
		vector<double> micro_dhpwl;
		vector<double> micro_doverlap;
		vector<double> micro_movement_window;
		vector<double> micro_dcongestion;
		vector<double> micro_overlap_history;
		vector<double> micro_overlap_x_history;
                vector<double> micro_overlap_x2_history;
		vector<double> micro_cdist_history;

		vector<double> accept_probabilities;

		void Logger::print_log(string logdir="./cache/") {
			std::ofstream clog(logdir+"cost_history");
			for(vector<double>::const_iterator i = cost_history.begin(); i!= cost_history.end(); ++i) {
				clog << *i << '\n';
			}
			clog.close();

			std::ofstream hlog(logdir+"hpwl_history");
			for(vector<double>::const_iterator i = hpwl_history.begin(); i!= hpwl_history.end(); ++i) {
				hlog << *i << '\n';
			}
			hlog.close();
			std::ofstream olog(logdir+"overlap_history");
			for(vector<double>::const_iterator i = overlap_history.begin(); i!= overlap_history.end(); ++i) {
				olog << *i << '\n';
			}
			std::ofstream xlog(logdir+"overlap_x_history");
			for(vector<double>::const_iterator i = overlap_x_history.begin(); i!= overlap_x_history.end(); ++i) {
				xlog << *i << '\n';
			}
			std::ofstream x2log(logdir+"overlap_x2_history");
			for(vector<double>::const_iterator i = overlap_x2_history.begin(); i!= overlap_x2_history.end(); ++i) {
				x2log << *i << '\n';
			}
			std::ofstream cdlog(logdir+"cdist_history");
			for(vector<double>::const_iterator i = cdist_history.begin(); i!= cdist_history.end(); ++i) {
				cdlog << *i << '\n';
			}
			olog.close();
			std::ofstream colog(logdir+"congestion_history");
			for(vector<double>::const_iterator i = congestion_history.begin(); i!= congestion_history.end(); ++i) {
				colog << *i << '\n';
			}
			colog.close();

			std::ofstream alog(logdir+"acceptance_ratio_history");
			for(vector<double>::const_iterator i = acceptance_ratio_history.begin(); i!= acceptance_ratio_history.end(); ++i) {
				alog << *i << '\n';
			}
			alog.close();
			std::ofstream tlog(logdir+"temperature_history");
			for(vector<double>::const_iterator i = temperature_history.begin(); i!= temperature_history.end(); ++i) {
				tlog << *i << '\n';
			}
			tlog.close();
			std::ofstream wlog(logdir+"move_window_coef_history");
			for(vector<double>::const_iterator i = move_window_coef_history.begin(); i!= move_window_coef_history.end(); ++i) {
				wlog << *i << '\n';
			}
			wlog.close();
			std::ofstream llog(logdir+"lambda_history");
	                for(vector<double>::const_iterator i = lambda_history.begin(); i!= lambda_history.end(); ++i) {
				llog << *i << '\n';
			}
			llog.close();
			std::ofstream lrlog(logdir+"lamrate_history");
                        for(vector<double>::const_iterator i = lamrate_history.begin(); i!= lamrate_history.end(); ++i) {
				lrlog << *i << '\n';
			}
			lrlog.close();

                        std::ofstream mdclog(logdir+"micro_dcost");
                        for(vector<double>::const_iterator i = micro_dcost.begin(); i!= micro_dcost.end(); ++i) {
				mdclog << *i << '\n';
			}
			mdclog.close();
                        std::ofstream mdhlog(logdir+"micro_dhpwl");
                        for(vector<double>::const_iterator i = micro_dhpwl.begin(); i!= micro_dhpwl.end(); ++i) {
				mdhlog << *i << '\n';
			}
			mdhlog.close();
                        std::ofstream mdolog(logdir+"micro_doverlap");
                        for(vector<double>::const_iterator i = micro_doverlap.begin(); i!= micro_doverlap.end(); ++i) {
				mdolog << *i << '\n';
			}
			mdolog.close();
                        std::ofstream mmlog(logdir+"micro_movement");
                        for(vector<double>::const_iterator i = micro_movement_window.begin(); i!= micro_movement_window.end(); ++i) {
				mmlog << *i << '\n';
			}
			mmlog.close();
                        std::ofstream mdcolog(logdir+"micro_dcongestion");
                        for(vector<double>::const_iterator i = micro_dcongestion.begin(); i!= micro_dcongestion.end(); ++i) {
				mdcolog << *i << '\n';
			}
			mdcolog.close();
			std::ofstream mxlog(logdir+"micro_overlap_x_history");
			for(vector<double>::const_iterator i = micro_overlap_x_history.begin(); i!= micro_overlap_x_history.end(); ++i) {
				mxlog << *i << '\n';
			}
			std::ofstream mx2log(logdir+"micro_overlap_x2_history");
			for(vector<double>::const_iterator i = micro_overlap_x2_history.begin(); i!= micro_overlap_x2_history.end(); ++i) {
				mx2log << *i << '\n';
			}
			std::ofstream mcdlog(logdir+"micro_cdist_history");
			for(vector<double>::const_iterator i = micro_cdist_history.begin(); i!= micro_cdist_history.end(); ++i) {
				mcdlog << *i << '\n';
			}

			std::ofstream mplog(logdir+"micro_accept_probs");
                        for(vector<double>::const_iterator i = accept_probabilities.begin(); i!= accept_probabilities.end(); ++i) {
				mplog << *i << '\n';
			}
			mplog.close();
		}
		void Logger::update_accept_probs(double prob) {
		    accept_probabilities.push_back(prob);
		}
		void Logger::update_micro_histories(double d_cost, double d_hpwl, double d_overlap, double d_congestion, double window) {
                    micro_dcost.push_back(d_cost);
		    micro_dhpwl.push_back(d_hpwl);
		    micro_doverlap.push_back(d_overlap);
		    micro_movement_window.push_back(window);
		    micro_dcongestion.push_back(d_congestion);
		}
		void Logger::update_micro_overlap_histories(double overlap_x, double overlap_x2, double cdist) {
                    micro_overlap_x_history.push_back(overlap_x);
		    micro_overlap_x2_history.push_back(overlap_x2);
		    micro_cdist_history.push_back(cdist);

		}
                void Logger::update_overlap_histories(double overlap_x, double overlap_x2, double cdist) {
                    overlap_x_history.push_back(overlap_x);
		    overlap_x2_history.push_back(overlap_x2);
		    cdist_history.push_back(cdist);

		}
		void Logger::update_cost_histories(double cost, double hpwl, double overlap, double congestion) {
			cost_history.push_back(cost);
			hpwl_history.push_back(hpwl);
			overlap_history.push_back(overlap);
			congestion_history.push_back(congestion);
		}
		void Logger::update_annealer_param_histories(double acceptance_ratio, double move_radius_coef, double temperature, double lambda, double lamrate) {
			acceptance_ratio_history.push_back(acceptance_ratio);
			move_window_coef_history.push_back(move_radius_coef);
			temperature_history.push_back(temperature);
			lambda_history.push_back(lambda);
			lamrate_history.push_back(lamrate);
		}
		void Logger::update_acceptance_ratio_history(double acceptance_ratio) {
			acceptance_ratio_history.push_back(acceptance_ratio);
		}
		void Logger::update_move_window_coef_history(double move_radius_coef) {
			move_window_coef_history.push_back(move_radius_coef);
		}
		void Logger::update_temperature_history(double temperature) {
			temperature_history.push_back(temperature);
		}
		void Logger::update_lambda_history(double lambda) {
			lambda_history.push_back(lambda);
		}
};

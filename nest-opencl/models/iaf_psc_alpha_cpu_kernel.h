#ifndef IAF_PSC_ALPHA_CPU_KERNEL_H
#define IAF_PSC_ALPHA_CPU_KERNEL_H

#include <vector>
using namespace std;

extern
void kernel_update(
  unsigned int tid,
  int num_nodes,
      long lag,
      int ring_buffer_size,
      ////
double *S__y3_,
double *P__Theta_,
double *V__P22_ex_,
double *V__P21_in_,
double *S__dI_ex_,
double *P__I_e_,
double *V__IPSCInitialValue_,
double *V__P31_ex_,
double *S__I_in_,
double *V__expm1_tau_m_,
int *S__r_,
double *S__I_ex_,
double *V__P21_ex_,
double *P__LowerBound_,
double *V__P22_in_,
double *V__weighted_spikes_ex_,
double *V__P11_in_,
double *V__weighted_spikes_in_,
double *V__P31_in_,
double *V__EPSCInitialValue_,
double *V__P32_ex_,
double *V__P11_ex_,
double *S__dI_in_,
double *P__V_reset_,
int *V__RefractoryCounts_,
double *V__P30_,
double *V__P32_in_,
double *S__y0_,
vector< double > &currents_,
vector< double > &ex_spikes_,
vector< double > &in_spikes_,
      ///
unsigned int *spike_count,
int time_index);

extern
void kernel_deliver_events(int num_nodes,
           int event_size,
           vector< int > &connections_ptr,
           vector< int > &connections,
           vector< double > &connections_weight,
           int num_spikes,
           int *spike_src,
           int *spike_multiplicity,
           long *spike_pos,
           vector< double > &ex_spikes_,
           vector< double > &in_spikes_,
           int time_index);

extern
void kernel_deliver_events_stdp_pl(
  unsigned int tid,
  int num_nodes, int batch_size,
             int ring_buffer_size,

            //  int *d_history_ptr,
            //  double *d_history_Kminus_,
            //  double *d_history_t_,
            //  int *d_history_access_counter_,
            //  double *d_Kminus_,
            //  double *d_tau_minus_inv_,

             vector< int > &d_history_ptr,
             vector< double > &d_history_Kminus_,
             vector< double > &d_history_t_,
             vector< int > &d_history_access_counter_,
             vector< double > &d_Kminus_,
             vector< double > &d_tau_minus_inv_,

             int *d_spike_tgid,
             double d_t_spike,
             double d_dendritic_delay,
             double *d_weight_,
             double *d_Kplus_,
             int *d_conn_type_,
             double *d_t_lastspike,
             long *d_pos,
             int *d_multiplicity,
             int event_multiplicity,
             int update_type,
             int d_conn_type,
             double cp_lambda,
             double cp_mu,
             double cp_alpha,
             double cp_tau_plus_inv,
             vector< double > &ex_spikes_,
             vector< double > &in_spikes_,
             int time_index);

#endif
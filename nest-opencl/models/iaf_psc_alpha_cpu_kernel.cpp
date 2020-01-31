#include "iaf_psc_alpha_cpu_kernel.h"

#include <cmath>
#include <iostream>
#include <assert.h>

using namespace std;

#define OPENCL_DBL_MEMCPY(dest, src, n, tid, nthreads)    \
    for (int i = 0; i < n; i++)         \
    {               \
      int idx = nthreads * i + tid;       \
      dest[idx] = src[idx];         \
    }

#define GSL_MAX_DBL(a,b) ((a) > (b) ? (a) : (b))

double ring_buffer_get_value(vector< double > &ring_buffer, int ring_buffer_size, int num_nodes, int tid, long lag, int time_index)
{
  int ind = num_nodes * ((time_index + lag) % ring_buffer_size) + tid;
  //if(ind >= ring_buffer_size) {
  //  cout << "ind " << ind << endl;
  //  cout << ring_buffer_size << " " << num_nodes << " " << tid << " " << lag << " " << time_index << endl;
    //exit(-1);
  //}
  //assert(ind < ring_buffer_size);
  double val = ring_buffer[ind];
  ring_buffer[ind] = 0.0;
  return val;
}

//void ring_buffer_add_value(double *ring_buffer, int ring_buffer_size, int num_nodes, int tgtid, long pos, double val, int time_index)
//{
//  int ind = num_nodes * ((time_index + pos) % ring_buffer_size) + tgtid;
//  ring_buffer[ind] += val;
//}

void ring_buffer_atomic_add_value(vector< double > &ring_buffer, int ring_buffer_size, int num_nodes, int tgtid, long pos, double val, int time_index)
{
  int ind = num_nodes * ((time_index + pos) % ring_buffer_size) + tgtid;
  //assert(ind < ring_buffer_size);
  ring_buffer[ind] += val;
  //AtomicAdd(ring_buffer + ind, val);
}

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
int time_index)
{
    //unsigned int tid = get_global_id(0);

  //if (tid >= num_nodes)
  //    return;

  {
    if ( S__r_[tid] == 0 )
    {
      // neuron not refractory
      S__y3_[tid] = V__P30_[tid] * ( S__y0_[tid] + P__I_e_[tid] ) + V__P31_ex_[tid] * S__dI_ex_[tid]
        + V__P32_ex_[tid] * S__I_ex_[tid] + V__P31_in_[tid] * S__dI_in_[tid] + V__P32_in_[tid] * S__I_in_[tid]
        + V__expm1_tau_m_[tid] * S__y3_[tid] + S__y3_[tid];

      // lower bound of membrane potential
      S__y3_[tid] = ( S__y3_[tid] < P__LowerBound_[tid] ? P__LowerBound_[tid] : S__y3_[tid] );
    }
    else
    {
      // neuron is absolute refractory
      --S__r_[tid];
    }

    // alpha shape EPSCs
    S__I_ex_[tid] = V__P21_ex_[tid] * S__dI_ex_[tid] + V__P22_ex_[tid] * S__I_ex_[tid];
    S__dI_ex_[tid] *= V__P11_ex_[tid];

    // Apply spikes delivered in this step; spikes arriving at T+1 have
    // an immediate effect on the state of the neuron
    V__weighted_spikes_ex_[tid] = ring_buffer_get_value(ex_spikes_, ring_buffer_size, num_nodes, tid, lag, time_index);
    S__dI_ex_[tid] += V__EPSCInitialValue_[tid] * V__weighted_spikes_ex_[tid];

    // alpha shape EPSCs
    S__I_in_[tid] = V__P21_in_[tid] * S__dI_in_[tid] + V__P22_in_[tid] * S__I_in_[tid];
    S__dI_in_[tid] *= V__P11_in_[tid];

    // Apply spikes delivered in this step; spikes arriving at T+1 have
    // an immediate effect on the state of the neuron
    V__weighted_spikes_in_[tid] = ring_buffer_get_value(in_spikes_, ring_buffer_size, num_nodes, tid, lag, time_index);
    S__dI_in_[tid] += V__IPSCInitialValue_[tid] * V__weighted_spikes_in_[tid];

    spike_count[tid] = 0;

    // threshold crossing
    if ( S__y3_[tid] >= P__Theta_[tid] )
    {
      S__r_[tid] = V__RefractoryCounts_[tid];
      S__y3_[tid] = P__V_reset_[tid];
      // A supra-threshold membrane potential should never be observable.
      // The reset at the time of threshold crossing enables accurate
      // integration independent of the computation step size, see [2,3] for
      // details.

      spike_count[tid]++;

      /*set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );
      SpikeEvent se;
      kernel().event_delivery_manager.send( *this, se, lag )*/
    }

    // set new input current
    S__y0_[tid] = ring_buffer_get_value(currents_, ring_buffer_size, num_nodes, tid, lag, time_index);

    // log state data
    /*LOG DATA*/
    //B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void insert_into_ring_buffer(int num_nodes,
                             int ring_buffer_size,
                             int tgtid,
                             long pos,
                             double weight,
                             int multiplicity,
                             vector< double > &ex_spikes_,
                             vector< double > &in_spikes_,
                             int time_index)
{
  const double s = weight * multiplicity;

  if ( weight > 0.0 )
  {
    ring_buffer_atomic_add_value(ex_spikes_, ring_buffer_size, num_nodes, tgtid, pos, s, time_index);
  }
  else
  {
    ring_buffer_atomic_add_value(in_spikes_, ring_buffer_size, num_nodes, tgtid, pos, s, time_index);
  }
}

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
           int time_index)
{
  //unsigned int tid = 0;//get_global_id(0);
  //unsigned int wid = tid / WARP_SIZE;
  //unsigned int lane = tid & (WARP_SIZE - 1);
  //unsigned int total_threads = 1;//get_global_size(0);
  //unsigned int total_warps = total_threads / WARP_SIZE;

  for (int src = 0; src < num_spikes; src ++)
  {
    int src_id = spike_src[src];
    int multiplicity = spike_multiplicity[src];
    long pos = spike_pos[src];
      
    int tgt_start = connections_ptr[src_id];
    int tgt_end = connections_ptr[src_id + 1];
      
    for (int i = tgt_start; i < tgt_end; i ++)
    {
      int tgt_id = connections[i];
      double weight = connections_weight[i];
      
      insert_into_ring_buffer(num_nodes,
            event_size,
            tgt_id,
            pos,
            weight,
            multiplicity,
            ex_spikes_,
            in_spikes_,
            time_index);
    }
  }
}

void get_history(double t1, double t2, int history_size, double* history_t_,
     int* history_access_counter_, int *start, int *finish)
{
  *finish = history_size;

  if (history_size == 0)
  {
    *start = *finish;
    return;
  }
  else
  {
    int runner = 0;
    while ((runner != history_size) && (history_t_[runner] <= t1))
    {
      ++runner;
    }
      *start = runner;
      while ((runner != history_size) && (history_t_[runner] <= t2))
  {
    history_access_counter_[runner]++;
    ++runner;
  }
      *finish = runner;
    }
}

double facilitate_(double w, double kplus, double lambda_, double mu_plus_, double Wmax_)
{
  double norm_w = ( w / Wmax_ )
    + ( lambda_ * pow( 1.0 - ( w / Wmax_ ), mu_plus_ ) * kplus );
  return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
}

double depress_( double w, double kminus, double lambda_, double alpha_, double mu_minus_, double Wmax_ )
{
  double norm_w = ( w / Wmax_ )
    - ( alpha_ * lambda_ * pow( w / Wmax_, mu_minus_ ) * kminus );
  return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
}

double facilitate_pl_(double w, double kplus, double cp_lambda_, double cp_mu_)
{
  return w + ( cp_lambda_ * pow( w, cp_mu_ ) * kplus );
}

double depress_pl_( double w, double kminus, double cp_lambda_, double cp_alpha_ )
{
  double new_w = w - ( cp_lambda_ * cp_alpha_ * w * kminus );
  return new_w > 0.0 ? new_w : 0.0;
}

double get_K_value(double t, int history_size,
       double K_minus_,
       double* history_Kminus_,
       double* history_t_,
       double tau_minus_inv_)
{
  if (history_size == 0)
    {
      return K_minus_;
    }
  int i = history_size - 1;
  while (i >= 0)
    {
      if (t > history_t_[i])
  {
    return (history_Kminus_[i] *
      exp((history_t_[i] - t) * tau_minus_inv_));
  }
      i--;
    }
  return 0;
}


void send_stdp(int num_nodes,
         int ring_buffer_size,
         int tgtid,
         double t_spike,
         double dendritic_delay,
         double *weight_,
         double *Kplus_,
         int conn_type_,
         double t_lastspike,
         long pos,
         //
         int history_size,
         double *history_Kminus_,
         double *history_t_,
         int *history_access_counter_,
         ///
         double K_minus_,
         double tau_minus_inv_,
         double cp_lambda,
         double cp_mu,
         double cp_alpha,
         double cp_tau_plus_inv,
         vector< double > &ex_spikes_,
         vector< double > &in_spikes_,
         int multiplicity,
         int time_index)
{
  if (conn_type_ == 2)
  {
    int start, finish;
    get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay, history_size,
    history_t_, history_access_counter_, &start, &finish);

    double minus_dt;
    while (start != finish)
    {
      minus_dt = t_lastspike - (history_t_[start] + dendritic_delay);
      start++;
      if (minus_dt == 0)
        {
          continue;
        }
      *weight_ = facilitate_pl_(*weight_, *Kplus_ * exp(minus_dt * cp_tau_plus_inv), cp_lambda, cp_mu);
    }

    double kminus = get_K_value(t_spike - dendritic_delay,  history_size,
        K_minus_,
        history_Kminus_,
        history_t_,
        tau_minus_inv_);
    *weight_ = depress_pl_( *weight_, kminus, cp_lambda, cp_alpha );

    insert_into_ring_buffer(num_nodes,
              ring_buffer_size,
              tgtid,
              pos,
              *weight_,
              multiplicity,
              ex_spikes_,
              in_spikes_,
              time_index);
        *Kplus_ = *Kplus_ * exp( ( t_lastspike - t_spike ) * cp_tau_plus_inv ) + 1.0;
  }
  else if (conn_type_ == 1)
  {
      insert_into_ring_buffer(num_nodes,
                  ring_buffer_size,
                  tgtid,
                  pos,
                  *weight_,
                  multiplicity,
                  ex_spikes_,
                  in_spikes_,
            time_index);
  }
}

void kernel_deliver_events_stdp_pl(
  unsigned int tid,
  int num_nodes, int batch_size,
             int ring_buffer_size,
             //
             vector< int > &d_history_ptr,
             vector< double > &d_history_Kminus_,
             vector< double > &d_history_t_,
             vector< int > &d_history_access_counter_,
             vector< double > &d_Kminus_,
             vector< double > &d_tau_minus_inv_,
             //
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
             int time_index)
{
  //unsigned int tid = get_global_id(0);

  //if (tid >= batch_size)
  //  return;
  
  int synapse_id = tid;
  int tgid = d_spike_tgid[synapse_id];
  double t_spike = d_t_spike;//[synapse_id];
  double dendritic_delay = d_dendritic_delay;//[synapse_id];
  double weight_ = d_weight_[synapse_id];
  double Kplus_ = d_Kplus_[synapse_id];
  int conn_type_ = d_conn_type == 1 ? 1 : d_conn_type_[synapse_id];
  double t_lastspike = d_t_lastspike[synapse_id];
  long pos = d_pos[synapse_id];
  int multiplicity = update_type == 2 ? d_multiplicity[synapse_id] : event_multiplicity;

  int history_index = d_history_ptr[tgid];
  int history_size = d_history_ptr[tgid + 1] - history_index;

  double *tgt_history_Kminus_ = (history_size == 0) ? nullptr : &d_history_Kminus_[0] + history_index;
  double *tgt_history_t_ = (history_size == 0) ? nullptr : &d_history_t_[0] + history_index;
  int *tgt_history_access_counter_ = (history_size == 0) ? nullptr : &d_history_access_counter_[0] + history_index;

  double K_minus_ = d_Kminus_[tgid];
  double tau_minus_inv_ = d_tau_minus_inv_[tgid];

  send_stdp(num_nodes,
      ring_buffer_size,
      tgid,
      t_spike,
      dendritic_delay,
      &weight_,
      &Kplus_,
      conn_type_,
      t_lastspike,
      pos,
      history_size,
      tgt_history_Kminus_,
      tgt_history_t_,
      tgt_history_access_counter_,
      K_minus_,
      tau_minus_inv_,
      cp_lambda,
      cp_mu,
      cp_alpha,
      cp_tau_plus_inv,
      ex_spikes_,
      in_spikes_,
      multiplicity,
      time_index);
  //t_lastspike = t_spike;

  if (conn_type_ == 2)
  {
    d_weight_[synapse_id] = weight_;
    d_Kplus_[synapse_id] = Kplus_;
  }
}
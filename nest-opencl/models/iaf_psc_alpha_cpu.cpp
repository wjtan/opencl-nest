#include "iaf_psc_alpha_cpu.h"
#include "iaf_psc_alpha_cpu_kernel.h"
#include "iaf_psc_alpha.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"
#include "vp_manager_impl.h"

// #define PROFILING_SIZE
#include "profile.h"

const size_t OFFSET = 0;

extern void writeNodes(nest::thread t, int count, double* nodes, size_t num_nodes);

nest::iaf_psc_alpha_cpu::iaf_psc_alpha_cpu()
  : is_data_ready( false )
  , graph_size( 0 )
  //, h_connections_ptr( NULL )
  //, h_connections( NULL )
  , event_size( 0 )
  , event_multiplicity( 1 )
  , ring_buffer_size( 0 )
  , is_initialized( false )
  , is_ring_buffer_ready( false )
  , is_gpu_initialized( false )
  , time_index (0)
 , h_S__y3_( NULL )
 , h_P__Theta_( NULL )
 , h_V__P22_ex_( NULL )
 , h_V__P21_in_( NULL )
 , h_S__dI_ex_( NULL )
 , h_P__I_e_( NULL )
 , h_V__IPSCInitialValue_( NULL )
 , h_V__P31_ex_( NULL )
 , h_S__I_in_( NULL )
 , h_V__expm1_tau_m_( NULL )
 , h_S__r_( NULL )
 , h_S__I_ex_( NULL )
 , h_V__P21_ex_( NULL )
 , h_P__LowerBound_( NULL )
 , h_V__P22_in_( NULL )
 , h_V__weighted_spikes_ex_( NULL )
 , h_V__P11_in_( NULL )
 , h_V__weighted_spikes_in_( NULL )
 , h_V__P31_in_( NULL )
 , h_V__EPSCInitialValue_( NULL )
 , h_V__P32_ex_( NULL )
 , h_V__P11_ex_( NULL )
 , h_S__dI_in_( NULL )
 , h_P__V_reset_( NULL )
 , h_V__RefractoryCounts_( NULL )
 , h_V__P30_( NULL )
 , h_V__P32_in_( NULL )
 , h_S__y0_( NULL )

//, h_currents_( NULL )
//, h_ex_spikes_( NULL )
//, h_in_spikes_( NULL )

//, h_currents__mark( NULL )
//, h_ex_spikes__mark( NULL )
//, h_in_spikes__mark( NULL )
//, h_currents__count( NULL )
//, h_ex_spikes__count( NULL )
//, h_in_spikes__count( NULL )
//, h_currents__index( NULL )
//, h_ex_spikes__index( NULL )
//, h_in_spikes__index( NULL )
  , h_spike_count (NULL)
  , isHistory (false)
//  , h_history_ptr (NULL)
//  , h_Kminus_ (NULL)
//  , h_tau_minus_inv_ (NULL)

{
}

nest::iaf_psc_alpha_cpu::~iaf_psc_alpha_cpu()
{
  //if (h_connections_ptr) delete[] h_connections_ptr;
  //if (h_connections) delete[] h_connections;
if (h_S__y3_) delete[] h_S__y3_;
if (h_P__Theta_) delete[] h_P__Theta_;
if (h_V__P22_ex_) delete[] h_V__P22_ex_;
if (h_V__P21_in_) delete[] h_V__P21_in_;
if (h_S__dI_ex_) delete[] h_S__dI_ex_;
if (h_P__I_e_) delete[] h_P__I_e_;
if (h_V__IPSCInitialValue_) delete[] h_V__IPSCInitialValue_;
if (h_V__P31_ex_) delete[] h_V__P31_ex_;
if (h_S__I_in_) delete[] h_S__I_in_;
if (h_V__expm1_tau_m_) delete[] h_V__expm1_tau_m_;
if (h_S__r_) delete[] h_S__r_;
if (h_S__I_ex_) delete[] h_S__I_ex_;
if (h_V__P21_ex_) delete[] h_V__P21_ex_;
if (h_P__LowerBound_) delete[] h_P__LowerBound_;
if (h_V__P22_in_) delete[] h_V__P22_in_;
if (h_V__weighted_spikes_ex_) delete[] h_V__weighted_spikes_ex_;
if (h_V__P11_in_) delete[] h_V__P11_in_;
if (h_V__weighted_spikes_in_) delete[] h_V__weighted_spikes_in_;
if (h_V__P31_in_) delete[] h_V__P31_in_;
if (h_V__EPSCInitialValue_) delete[] h_V__EPSCInitialValue_;
if (h_V__P32_ex_) delete[] h_V__P32_ex_;
if (h_V__P11_ex_) delete[] h_V__P11_ex_;
if (h_S__dI_in_) delete[] h_S__dI_in_;
if (h_P__V_reset_) delete[] h_P__V_reset_;
if (h_V__RefractoryCounts_) delete[] h_V__RefractoryCounts_;
if (h_V__P30_) delete[] h_V__P30_;
if (h_V__P32_in_) delete[] h_V__P32_in_;
if (h_S__y0_) delete[] h_S__y0_;
  if (h_spike_count) delete[] h_spike_count;
//if (h_currents_) delete[] h_currents_;
//if (h_ex_spikes_) delete[] h_ex_spikes_;
//if (h_in_spikes_) delete[] h_in_spikes_;


//if (h_currents__mark) delete[] h_currents__mark;
//if (h_ex_spikes__mark) delete[] h_ex_spikes__mark;
//if (h_in_spikes__mark) delete[] h_in_spikes__mark;
//if (h_currents__count) delete[] h_currents__count;
//if (h_ex_spikes__count) delete[] h_ex_spikes__count;
//if (h_in_spikes__count) delete[] h_in_spikes__count;
//if (h_currents__index) delete[] h_currents__index;
//if (h_ex_spikes__index) delete[] h_ex_spikes__index;
//if (h_in_spikes__index) delete[] h_in_spikes__index;
}

void
nest::iaf_psc_alpha_cpu::initialize_gpu()
{
  std::cout << "iaf_psc_alpha_cpu::initialize_gpu" << std::endl;
  if (not is_gpu_initialized)
    {
      //int thrd = kernel().vp_manager.get_thread_id();
      //if (thrd == 0)
      //{
      //  if (initialize_opencl_context())
      //    return;
      //}

      this->total_num_nodes = kernel().node_manager.size();
      connections.resize(total_num_nodes);
      is_gpu_initialized = true;
    }
}

void
nest::iaf_psc_alpha_cpu::initialize_nodes(const std::vector< Node* > &nodes) {
  const int t = kernel().vp_manager.get_thread_id();

  const Token model = kernel().model_manager.get_modeldict()->lookup( "iaf_psc_alpha" );
  const index model_id = static_cast< index >( model );

  for (std::vector< Node* >::const_iterator node = nodes.begin(); node != nodes.end(); ++node ){
    if ((*node)->get_model_id() == model_id) {
      actualNodes.push_back(*node);
    } else {
      otherNodes.push_back(*node);
    }
  }
  this->num_local_nodes = actualNodes.size();

  // //FIXME
  // total_num_nodes = actualNodes.size();
  // connections.resize(total_num_nodes);

  // BuildGraphEvent e;
  // for (size_t nodeId = 0; nodeId < actualNodes.size(); nodeId++) {
  //   e.set_sender_gid( nodeId );
  //   kernel().connection_manager.send( t, nodeId, e );
  // }

  std::cout << "Other Nodes " << otherNodes.size() << std::endl;
}

void
nest::iaf_psc_alpha_cpu::mass_update(const std::vector< Node* > &nodes, Time const& origin,
                    const long from,
                    const long to )
{
  mass_update_(nodes, origin, from, to, false );
}

bool
nest::iaf_psc_alpha_cpu::mass_wfr_update(const std::vector< Node* > &nodes, Time const& origin,
                        const long from,
                        const long to )
{
  return mass_update_(nodes, origin, from, to, true );
}

int gpuNodeCount = 0;

// TODO: the for loops will later be kernel calls
bool
nest::iaf_psc_alpha_cpu::mass_update_( const std::vector<Node *> &nodes,
                      Time const& origin,
                      const long from,
                      const long to,
                      const bool called_from_wfr_update ) // TODO: don't know yet whether we need to cover both cases here
{
  // TODO: do AoS for now, SoA will come later on

  PROFILING_INIT();

  //#ifdef PROFILING
  int thrd_id = kernel().vp_manager.get_thread_id();
  //#endif  

  //int count_spike = 0;
  
  initialize_device();
  
  // TODO: for now, I assume that these are the same for all nodes

  if (not is_data_ready)
  {
    cout << "Copy Data to Device" << endl;
    //copy_data_to_device(nodes);
    copy_data_to_device(actualNodes);
    is_data_ready = true;
  }
  // if (from < to)
  //   prepare_copy_to_device(nodes, called_from_wfr_update, from);
  
  //set_kernel_args(this->gpu_kernel, this->num_local_nodes);

  for ( long lag = from; lag < to; ++lag)
  {
  //     for ( std::vector<Node*>::iterator nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++ )
  //     {
  //       nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*nodeIt;
  //       node->pre_gsl(origin, lag);
  //     }

    //set_lag_args(this->gpu_kernel, lag);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Start of GPU section

    //PROFILING_START();

    //fill_buffer_zero_uint(&gpu_context, d_spike_count, this->num_local_nodes*sizeof(unsigned int));
    for(size_t i = 0; i < this->num_local_nodes; i++){
      h_spike_count[i] = 0;
    }
    //synchronize();
    //this->command_queue.flush();

    //PROFILING_END("HtD");

    PROFILING_START();
    //execute_kernel(this->gpu_kernel, &gpu_context, this->num_local_nodes);
    

    for(size_t tid = 0; tid < this->num_local_nodes - OFFSET; tid++) {

kernel_update(
  tid,
  this->num_local_nodes - OFFSET,
      lag,
  this->event_size,
      ////
h_S__y3_,
h_P__Theta_,
h_V__P22_ex_,
h_V__P21_in_,
h_S__dI_ex_,
h_P__I_e_,
h_V__IPSCInitialValue_,
h_V__P31_ex_,
h_S__I_in_,
h_V__expm1_tau_m_,
h_S__r_,
h_S__I_ex_,
h_V__P21_ex_,
h_P__LowerBound_,
h_V__P22_in_,
h_V__weighted_spikes_ex_,
h_V__P11_in_,
h_V__weighted_spikes_in_,
h_V__P31_in_,
h_V__EPSCInitialValue_,
h_V__P32_ex_,
h_V__P11_ex_,
h_S__dI_in_,
h_P__V_reset_,
h_V__RefractoryCounts_,
h_V__P30_,
h_V__P32_in_,
h_S__y0_,
h_currents_,
h_ex_spikes_,
h_in_spikes_,
      ///
h_spike_count,
time_index);
}

    // if (lag + 1 < to)
    //     prepare_copy_to_device(nodes, called_from_wfr_update, lag + 1);

    //synchronize();
    //this->command_queue.flush();
    PROFILING_END("Execute kernel");

    //writeNodes(thrd_id, gpuNodeCount++, h_S__y3_, this->num_local_nodes);

    //PROFILING_START();
    //copy_data_from_device(nodes, false);
    //synchronize();
    //PROFILING_END("DtH");

    // End of GPU section
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    PROFILING_START();

    int count_spike = 0;
    int node_id = 0;
    //for ( std::vector<Node*>::const_iterator nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++, node_id++ )
    for ( std::vector<Node*>::const_iterator nodeIt = actualNodes.begin(); nodeIt != actualNodes.end(); nodeIt++, node_id++ )
    {
      //if (node_id >= this->num_local_nodes - OFFSET)
      //{
      //  break;
      //}

      nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*nodeIt;
      //node->post_gsl(origin, lag);

      if (h_spike_count[node_id] > 0)
      //for (size_t i = 0; i < h_spike_count[node_id]; i++)
      {
        node->set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );
        SpikeEvent se;
        kernel().event_delivery_manager.send( *node, se, lag );
        count_spike += h_spike_count[node_id];
        
        //kernel().simulation_manager.incSpikes();
      }
    }
    cout << lag << ") SpikeCount " << count_spike << endl;
    PROFILING_END("Spike send");
  }

  // PROFILING_START();
  // for (int nodeid = this->num_local_nodes - OFFSET; nodeid < nodes.size(); nodeid++)
  //   if (called_from_wfr_update) {
  //     nodes[nodeid]->wfr_update(origin, from, to);
  //   } else {
  //     nodes[nodeid]->update(origin, from, to);
  //   }
  // PROFILING_END("Node Update");

  if (called_from_wfr_update) {
    for( std::vector<Node*>::const_iterator node = otherNodes.begin(); node != otherNodes.end(); node++){
      (*node)->wfr_update(origin, from, to);
    }
  } else {
    for( std::vector<Node*>::const_iterator node = otherNodes.begin(); node != otherNodes.end(); node++){
      (*node)->update(origin, from, to);
    }
  }

  //PROFILING_START();
  //copy_data_from_device(nodes, true);
  //synchronize();
  //PROFILING_END("Last DtH");

  return true;
}

// int
// nest::iaf_psc_alpha_cpu::initialize_opencl_context()
// {
//   int thrd_id = kernel().vp_manager.get_thread_id();
//   printf("[%d] Initialize OpenCL Context\n", thrd_id);

//   const int n_gpus = kernel().vp_manager.get_num_gpus();
//   try
//     {
//       // Get list of OpenCL platforms.
//       std::vector<cl::Platform> list_platform;
//       cl::Platform::get(&list_platform);

//       if (list_platform.empty()) {
//         std::cerr << "OpenCL platforms not found." << std::endl;
//         return 1;
//       }

//       //gpu_context.platform = list_platform[0];
//       if (list_platform.size() == 1) {
//         gpu_context.platform = list_platform[0];
//       } else {
//         gpu_context.platform = list_platform[1];
//       }

//       std::cout << "Select Platform: " << gpu_context.platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

//       gpu_context.platform.getDevices(CL_DEVICE_TYPE_GPU, &gpu_context.list_device);
//       //gpu_context.platform.getDevices(CL_DEVICE_TYPE_CPU, &gpu_context.list_device);

//       context = cl::Context(gpu_context.list_device);

//       string source = "iaf_psc_alpha";
//       string kernel_path = "kernel/";
//       string src_file = kernel_path + source + ".cl";
//       FILE *fs = fopen(src_file.c_str(), "r");
//       if (!fs)
//       {
//         std::cerr << "Failed to load kernel file." << std::endl;
//         return 1;
//       }

//       char *source_str = (char *)malloc(MAX_SOURCE_SIZE);
//       size_t source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fs);
//       fclose(fs);

//       // Compile OpenCL program for found device.
//       program = cl::Program(context,
//                 cl::Program::Sources(1, std::make_pair(source_str, source_size)));

//       try {
//         program.build(gpu_context.list_device, "-cl-nv-maxrregcount=200 -cl-nv-verbose");
//         //program.build(gpu_context.list_device, "-cl-fast-relaxed-math");
//         //program.build(gpu_context.list_device);
//       } catch (const cl::Error&) {
//         std::cerr
//           << "OpenCL compilation error" << std::endl
//           << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(gpu_context.list_device[0])
//           << std::endl;
//         return 1;
//       }

//       const int n_devices = gpu_context.list_device.size();
      
//       if (n_devices < n_gpus) {
//         kernel().vp_manager.set_num_gpus(n_devices);
//         std::cout << "Set num of GPUS=" << n_devices << " (" << n_gpus << ")" << std::endl;
//       }

//     }
//     catch (const cl::Error &err) {
//       std::cerr
//           << "OpenCL error: "
//           << err.what() << "(" << err.err() << ")"
//           << std::endl;
//       return 1;
//     }

//   return 0;
// }

// int
// nest::iaf_psc_alpha_cpu::initialize_command_queue()
// {
//   const int num_gpus = kernel().vp_manager.get_num_gpus();
//   const int thrd_id = kernel().vp_manager.get_thread_id();
//   const int vp_id = kernel().vp_manager.thread_to_vp(thrd_id);

//   try
//   {
//     //printf("[%d] Num of Devices: %d\n", thrd_id, num_devices);

//     this->command_queue = cl::CommandQueue(context, gpu_context.list_device[vp_id % num_gpus]);
//   }
//   catch (const cl::Error &err) {
//     std::cerr
//         << "OpenCL error: "
//         << err.what() << "(" << err.err() << ")"
//         << std::endl;
//     return 1;
//   }

//   return 0;
// }

void
nest::iaf_psc_alpha_cpu::initialize_device()
{
  if (!is_initialized)
    {
      cout << "initialize_device" << endl;
      int thrd_id = kernel().vp_manager.get_thread_id();
      int len = this->num_local_nodes;
      //setlocale(LC_NUMERIC, "");
      printf("[%d] Num of nodes: %ld\n", thrd_id, this->num_local_nodes);
      //printf("[%d] Double: %'ld Int: %'ld\n", thrd_id, 26 * len * sizeof(double), 3 * len * sizeof(int));
      
      h_S__y3_ = new double[len];
      h_P__Theta_ = new double[len];
      h_V__P22_ex_ = new double[len];
      h_V__P21_in_ = new double[len];
      h_S__dI_ex_ = new double[len];
      h_P__I_e_ = new double[len];
      h_V__IPSCInitialValue_ = new double[len];
      h_V__P31_ex_ = new double[len];
      h_S__I_in_ = new double[len];
      h_V__expm1_tau_m_ = new double[len];
      h_S__r_ = new int[len];
      h_S__I_ex_ = new double[len];
      h_V__P21_ex_ = new double[len];
      h_P__LowerBound_ = new double[len];
      h_V__P22_in_ = new double[len];
      h_V__weighted_spikes_ex_ = new double[len];
      h_V__P11_in_ = new double[len];
      h_V__weighted_spikes_in_ = new double[len];
      h_V__P31_in_ = new double[len];
      h_V__EPSCInitialValue_ = new double[len];
      h_V__P32_ex_ = new double[len];
      h_V__P11_ex_ = new double[len];
      h_S__dI_in_ = new double[len];
      h_P__V_reset_ = new double[len];
      h_V__RefractoryCounts_ = new int[len];
      h_V__P30_ = new double[len];
      h_V__P32_in_ = new double[len];
      h_S__y0_ = new double[len];
      h_spike_count = new unsigned int[len];

      // create(&gpu_context, &S__y3_, len*sizeof(double));
      // create(&gpu_context, &P__Theta_, len*sizeof(double));
      // create(&gpu_context, &V__P22_ex_, len*sizeof(double));
      // create(&gpu_context, &V__P21_in_, len*sizeof(double));
      // create(&gpu_context, &S__dI_ex_, len*sizeof(double));
      // create(&gpu_context, &P__I_e_, len*sizeof(double));
      // create(&gpu_context, &V__IPSCInitialValue_, len*sizeof(double));
      // create(&gpu_context, &V__P31_ex_, len*sizeof(double));
      // create(&gpu_context, &S__I_in_, len*sizeof(double));
      // create(&gpu_context, &V__expm1_tau_m_, len*sizeof(double));
      // create(&gpu_context, &S__r_, len*sizeof(int));
      // create(&gpu_context, &S__I_ex_, len*sizeof(double));
      // create(&gpu_context, &V__P21_ex_, len*sizeof(double));
      // create(&gpu_context, &P__LowerBound_, len*sizeof(double));
      // create(&gpu_context, &V__P22_in_, len*sizeof(double));
      // create(&gpu_context, &V__weighted_spikes_ex_, len*sizeof(double));
      // create(&gpu_context, &V__P11_in_, len*sizeof(double));
      // create(&gpu_context, &V__weighted_spikes_in_, len*sizeof(double));
      // create(&gpu_context, &V__P31_in_, len*sizeof(double));
      // create(&gpu_context, &V__EPSCInitialValue_, len*sizeof(double));
      // create(&gpu_context, &V__P32_ex_, len*sizeof(double));
      // create(&gpu_context, &V__P11_ex_, len*sizeof(double));
      // create(&gpu_context, &S__dI_in_, len*sizeof(double));
      // create(&gpu_context, &P__V_reset_, len*sizeof(double));
      // create(&gpu_context, &V__RefractoryCounts_, len*sizeof(int));
      // create(&gpu_context, &V__P30_, len*sizeof(double));
      // create(&gpu_context, &V__P32_in_, len*sizeof(double));
      // create(&gpu_context, &S__y0_, len*sizeof(double));
      // create(&gpu_context, &d_spike_count, len*sizeof(unsigned int));

      is_initialized = true;
    }
}

#define UPLOAD_2D_DATA(buf, src)        \
  if (node->src != NULL)            \
    for (int j = 0; j < dimension; j++)        \
      buf[j*num_nodes + i] = node->src[j];                

#define UPLOAD_1D_DATA(buf, src)        \
  buf[i] = node->src;                            

#define UPLOAD_1D_DATA_VALUE(buf, val)        \
  buf[i] = val;                            

//#define FINISH_2D_UPLOAD(dst, buf, dim, type)                \
//  upload(&gpu_context, (void *)buf, dst, dim*num_nodes*sizeof(type));

//#define FINISH_1D_UPLOAD(dst, buf, type)                \
//  upload(&gpu_context, (void *)buf, dst, num_nodes*sizeof(type));

// void
// nest::iaf_psc_alpha_cpu::prepare_copy_to_device(std::vector< Node* > &nodes, bool called_from_wfr_update, long lag_)
// {
  // int num_nodes = nodes.size();

  // int wfr_interpolation_order = kernel().simulation_manager.get_wfr_interpolation_order();

  // //int len3 = /*buffer_size*/ 4 * num_nodes;
    
  // std::vector<Node *>::iterator nodeIt = nodes.begin();

  // for (int i = 0; nodeIt != nodes.end(); nodeIt++, i++ )
  //   {
  //     nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*nodeIt;
      
  //     double *tmp = h_B_interpolation_coefficients;// + len3;

  //     switch (wfr_interpolation_order)
  //     {
  //     case 0:
  //       tmp[num_nodes*0 + i] = node->B_.interpolation_coefficients[ lag_ ];
  //       break;
  //     case 1:
  //       tmp[num_nodes*0 + i] = node->B_.interpolation_coefficients[ lag_ * 2 + 0 ];
  //       tmp[num_nodes*1 + i] = node->B_.interpolation_coefficients[ lag_ * 2 + 1 ];
  //       break;
  //     case 3:
  //       tmp[num_nodes*0 + i] = node->B_.interpolation_coefficients[ lag_ * 4 + 0 ];
  //       tmp[num_nodes*1 + i] = node->B_.interpolation_coefficients[ lag_ * 4 + 1 ];
  //       tmp[num_nodes*2 + i] = node->B_.interpolation_coefficients[ lag_ * 4 + 2 ];
  //       tmp[num_nodes*3 + i] = node->B_.interpolation_coefficients[ lag_ * 4 + 3 ];
  //       break;
  //     default:
  //       break;
  //     }

  //     if (called_from_wfr_update)
  //     {
  //       double val = node->B_.spike_exc_.get_value_wfr_update( lag_ ) * node->V_.PSCurrInit_E_;
  //       UPLOAD_1D_DATA_VALUE(h_B_spike_exc_, val);
  //       val = node->B_.spike_inh_.get_value_wfr_update( lag_ ) * node->V_.PSCurrInit_I_;
  //       UPLOAD_1D_DATA_VALUE(h_B_spike_inh_, val);
  //     }
  //   }
// }

void
nest::iaf_psc_alpha_cpu::copy_data_to_device(const std::vector< Node* > &nodes)
{
  int num_nodes = nodes.size();

  std::vector<Node *>::const_iterator nodeIt = nodes.begin();

  for (int i = 0; nodeIt != nodes.end(); nodeIt++, i++ )
  {
    nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*nodeIt;

UPLOAD_1D_DATA(h_S__y3_, S_.y3_);
UPLOAD_1D_DATA(h_P__Theta_, P_.Theta_);
UPLOAD_1D_DATA(h_V__P22_ex_, V_.P22_ex_);
UPLOAD_1D_DATA(h_V__P21_in_, V_.P21_in_);
UPLOAD_1D_DATA(h_S__dI_ex_, S_.dI_ex_);
UPLOAD_1D_DATA(h_P__I_e_, P_.I_e_);
UPLOAD_1D_DATA(h_V__IPSCInitialValue_, V_.IPSCInitialValue_);
UPLOAD_1D_DATA(h_V__P31_ex_, V_.P31_ex_);
UPLOAD_1D_DATA(h_S__I_in_, S_.I_in_);
UPLOAD_1D_DATA(h_V__expm1_tau_m_, V_.expm1_tau_m_);
UPLOAD_1D_DATA(h_S__r_, S_.r_);
UPLOAD_1D_DATA(h_S__I_ex_, S_.I_ex_);
UPLOAD_1D_DATA(h_V__P21_ex_, V_.P21_ex_);
UPLOAD_1D_DATA(h_P__LowerBound_, P_.LowerBound_);
UPLOAD_1D_DATA(h_V__P22_in_, V_.P22_in_);
UPLOAD_1D_DATA(h_V__weighted_spikes_ex_, V_.weighted_spikes_ex_);
UPLOAD_1D_DATA(h_V__P11_in_, V_.P11_in_);
UPLOAD_1D_DATA(h_V__weighted_spikes_in_, V_.weighted_spikes_in_);
UPLOAD_1D_DATA(h_V__P31_in_, V_.P31_in_);
UPLOAD_1D_DATA(h_V__EPSCInitialValue_, V_.EPSCInitialValue_);
UPLOAD_1D_DATA(h_V__P32_ex_, V_.P32_ex_);
UPLOAD_1D_DATA(h_V__P11_ex_, V_.P11_ex_);
UPLOAD_1D_DATA(h_S__dI_in_, S_.dI_in_);
UPLOAD_1D_DATA(h_P__V_reset_, P_.V_reset_);
UPLOAD_1D_DATA(h_V__RefractoryCounts_, V_.RefractoryCounts_);
UPLOAD_1D_DATA(h_V__P30_, V_.P30_);
UPLOAD_1D_DATA(h_V__P32_in_, V_.P32_in_);
UPLOAD_1D_DATA(h_S__y0_, S_.y0_);


      /* DEVICE OUTPUT VAR UPLOAD */

    }
  

// FINISH_1D_UPLOAD(S__y3_, h_S__y3_, double);
// FINISH_1D_UPLOAD(P__Theta_, h_P__Theta_, double);
// FINISH_1D_UPLOAD(V__P22_ex_, h_V__P22_ex_, double);
// FINISH_1D_UPLOAD(V__P21_in_, h_V__P21_in_, double);
// FINISH_1D_UPLOAD(S__dI_ex_, h_S__dI_ex_, double);
// FINISH_1D_UPLOAD(P__I_e_, h_P__I_e_, double);
// FINISH_1D_UPLOAD(V__IPSCInitialValue_, h_V__IPSCInitialValue_, double);
// FINISH_1D_UPLOAD(V__P31_ex_, h_V__P31_ex_, double);
// FINISH_1D_UPLOAD(S__I_in_, h_S__I_in_, double);
// FINISH_1D_UPLOAD(V__expm1_tau_m_, h_V__expm1_tau_m_, double);
// FINISH_1D_UPLOAD(S__r_, h_S__r_, int);
// FINISH_1D_UPLOAD(S__I_ex_, h_S__I_ex_, double);
// FINISH_1D_UPLOAD(V__P21_ex_, h_V__P21_ex_, double);
// FINISH_1D_UPLOAD(P__LowerBound_, h_P__LowerBound_, double);
// FINISH_1D_UPLOAD(V__P22_in_, h_V__P22_in_, double);
// FINISH_1D_UPLOAD(V__weighted_spikes_ex_, h_V__weighted_spikes_ex_, double);
// FINISH_1D_UPLOAD(V__P11_in_, h_V__P11_in_, double);
// FINISH_1D_UPLOAD(V__weighted_spikes_in_, h_V__weighted_spikes_in_, double);
// FINISH_1D_UPLOAD(V__P31_in_, h_V__P31_in_, double);
// FINISH_1D_UPLOAD(V__EPSCInitialValue_, h_V__EPSCInitialValue_, double);
// FINISH_1D_UPLOAD(V__P32_ex_, h_V__P32_ex_, double);
// FINISH_1D_UPLOAD(V__P11_ex_, h_V__P11_ex_, double);
// FINISH_1D_UPLOAD(S__dI_in_, h_S__dI_in_, double);
// FINISH_1D_UPLOAD(P__V_reset_, h_P__V_reset_, double);
// FINISH_1D_UPLOAD(V__RefractoryCounts_, h_V__RefractoryCounts_, int);
// FINISH_1D_UPLOAD(V__P30_, h_V__P30_, double);
// FINISH_1D_UPLOAD(V__P32_in_, h_V__P32_in_, double);
// FINISH_1D_UPLOAD(S__y0_, h_S__y0_, double);

    /* DEVICE OUTPUT VAR FINISH UPLOAD */

  
}

// #define START_2D_DOWNLOAD(src, buf, dim, type)                \
//   download(&gpu_context, src, (void *)buf, dim*num_nodes*sizeof(type));

#define DOWNLOAD_2D_DATA(dst, buf)        \
  for (int j = 0; j < dimension; j++)        \
    node->dst[j] = buf[j*num_nodes + i];

// #define START_1D_DOWNLOAD(src, buf, type)                \
//   download(&gpu_context, src, (void *)buf, num_nodes*sizeof(type));

// #define START_1D_DOWNLOAD_OFF(src, buf, type, off)                \
//   download(&gpu_context, src, (void *)buf, num_nodes*sizeof(type), off);

#define DOWNLOAD_1D_DATA(dst, buf)        \
  node->dst = buf[i];                            

//void
//nest::iaf_psc_alpha_cpu::copy_data_from_device(const std::vector< Node* > &nodes, bool last_copy)
//{
//  int num_nodes = nodes.size();

//  std::vector<Node *>::const_iterator nodeIt = nodes.begin();

  // if (last_copy)
  //   {
  //     START_1D_DOWNLOAD(B_IntegrationStep_, h_B_IntegrationStep_, double);
  //     /* DEVICE OUTPUT VAR DOWNLOAD */
  //   }
    
  //START_1D_DOWNLOAD(d_spike_count, h_spike_count, unsigned int);
  
  // for (int i = 0; nodeIt != nodes.end(); nodeIt++, i++)
  //   {
  //     nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*nodeIt;

  //     if (last_copy)
  //     {
  //       DOWNLOAD_1D_DATA(B_.IntegrationStep_, h_B_IntegrationStep_);
  //     }

  //     /* DEVICE OUTPUT VAR FINISH DOWNLOAD */
     // }
//}

// void
// nest::iaf_psc_alpha_cpu::create(clContext_ *clCxt, cl::Buffer *mem, int len)
// {
//   cl_int ret;
//   //printf("len %d\n", len);
//   if (len == 0)
//     return;
//   *mem = cl::Buffer(context, CL_MEM_READ_WRITE, len, NULL, &ret);
//   if(ret != CL_SUCCESS){
//     printf("Failed to create buffer on GPU. %d\n", ret);
//     return ;
//   }
// }

// void
// nest::iaf_psc_alpha_cpu::upload(clContext_ *clCxt, void *data, cl::Buffer &gdata, int datalen)
// {
//   //write data to buffer
//   //printf("write datalen %d\n", datalen);
//   if (datalen == 0)
//     return;

//   try {
//     cl_int ret = this->command_queue.enqueueWriteBuffer(gdata, CL_FALSE, 0, datalen, (void *)data, NULL, NULL);
//     if(ret != CL_SUCCESS){
//       std::cerr << "OpenCL clEnqueueWriteBuffer failed." << ret << std::endl;
//       return ;
//     }
//   } catch (cl::Error err) {
//     std::cerr << "OpenCL enqueueWriteBuffer " << err.what() << "(" << err.err() << ")" << std::endl;
//     return;
//    }
// }

// void
// nest::iaf_psc_alpha_cpu::download(clContext_ *clCxt, cl::Buffer &gdata,void *data,int data_len, int offset)
// {
//   cl_int ret;
//   //printf("read datalen %d\n", data_len);
//   if (data_len == 0)
//     return;
//   ret = this->command_queue.enqueueReadBuffer(gdata, CL_FALSE, offset, data_len, (void *)data, NULL,NULL);
//   if(ret != CL_SUCCESS){
//     printf("clEnqueueReadBuffer failed. error code: %d\n", ret);
//     return ;
//   }
// }

// inline void
// nest::iaf_psc_alpha_cpu::synchronize()
// {
//   this->command_queue.finish();
// }

// int
// nest::iaf_psc_alpha_cpu::getKernel(string kernelName, string deliver_kernel_name, string static_deliver_kernel_name, clContext_ *clCxt)
// {
//   this->gpu_kernel = new cl::Kernel(program, kernelName.c_str());

//   this->deliver_kernel = new cl::Kernel(program, deliver_kernel_name.c_str());

//   this->static_deliver_kernel = new cl::Kernel(program, static_deliver_kernel_name.c_str());
  
//   return 0;
// }


// #define set_kernel(value)                        \
//   ret = kernel->setArg(arg_idx, value);                    \
//   if (ret != CL_SUCCESS) {                        \
//     printf("Failed to set arg %d, error code %d\n", arg_idx, ret);    \
//     return 1;                                \
//   }                                    \
//   arg_idx++;

// #define set_kernel_prim(value, type)                    \
//   ret = kernel->setArg(arg_idx, static_cast<type>(value));        \
//   if (ret != CL_SUCCESS) {                        \
//     printf("Failed to set arg %d, error code %d\n", arg_idx, ret);    \
//     return 1;                                \
//   }                                    \
//   arg_idx++;

// int
// nest::iaf_psc_alpha_cpu::set_lag_args(cl::Kernel *kernel, long lag)
// {
//   cl_int ret;
//   ret = kernel->setArg(1, static_cast<cl_long>(lag));
//   if (ret != CL_SUCCESS) {                        
//     printf("Failed to set lag arg, error code %dn", ret);    
//     return 1;                                
//   }                                    
//   return 0;
// }

// int
// nest::iaf_psc_alpha_cpu::set_kernel_args(cl::Kernel *kernel, int num_nodes)
// {
//   cl_int ret;
//   int arg_idx = 0;

//   set_kernel_prim(num_nodes, cl_int);

//   //set_kernel_prim(lag, cl_long);
//   arg_idx++;

//   set_kernel_prim(event_size, cl_int);

//   set_kernel(S__y3_);
//   set_kernel(P__Theta_);
//   set_kernel(V__P22_ex_);
//   set_kernel(V__P21_in_);
//   set_kernel(S__dI_ex_);
//   set_kernel(P__I_e_);
//   set_kernel(V__IPSCInitialValue_);
//   set_kernel(V__P31_ex_);
//   set_kernel(S__I_in_);
//   set_kernel(V__expm1_tau_m_);
//   set_kernel(S__r_);
//   set_kernel(S__I_ex_);
//   set_kernel(V__P21_ex_);
//   set_kernel(P__LowerBound_);
//   set_kernel(V__P22_in_);
//   set_kernel(V__weighted_spikes_ex_);
//   set_kernel(V__P11_in_);
//   set_kernel(V__weighted_spikes_in_);
//   set_kernel(V__P31_in_);
//   set_kernel(V__EPSCInitialValue_);
//   set_kernel(V__P32_ex_);
//   set_kernel(V__P11_ex_);
//   set_kernel(S__dI_in_);
//   set_kernel(P__V_reset_);
//   set_kernel(V__RefractoryCounts_);
//   set_kernel(V__P30_);
//   set_kernel(V__P32_in_);
//   set_kernel(S__y0_);
//   set_kernel(d_currents_);
//   set_kernel(d_ex_spikes_);
//   set_kernel(d_in_spikes_);
//   set_kernel(d_spike_count);
  
//   set_kernel_prim(time_index, cl_int);

//   return 0;
// }

// int
// nest::iaf_psc_alpha_cpu::set_deliver_kernel_args(cl::Kernel *kernel, int num_nodes, int batch_size)
// {
//   cl_int ret;
//   int arg_idx = 0;

//   set_kernel_prim(num_nodes, cl_int);
//   set_kernel_prim(batch_size, cl_int);
//   set_kernel_prim(event_size, cl_int);
//   set_kernel(d_history_ptr);
//   set_kernel(d_history_Kminus_);
//   set_kernel(d_history_t_);
//   set_kernel(d_history_access_counter_);
//   set_kernel(d_Kminus_);
//   set_kernel(d_tau_minus_inv_);
//   set_kernel(d_spike_tgid);
//   set_kernel_prim(event_t_spike, cl_double);
//   set_kernel_prim(event_dendritic_delay, cl_double);
//   set_kernel(d_weight_);
//   set_kernel(d_Kplus_);
//   set_kernel(d_conn_type_);
//   set_kernel(d_t_lastspike);
//   set_kernel(d_pos);
//   set_kernel(d_multiplicity);
//   set_kernel_prim(event_multiplicity, cl_int);
//   set_kernel_prim(this->update_type, cl_int);
//   set_kernel_prim(conn_type, cl_int);
//   set_kernel_prim(cp_lambda_, cl_double);
//   set_kernel_prim(cp_mu_, cl_double);
//   set_kernel_prim(cp_alpha_, cl_double);
//   set_kernel_prim(cp_tau_plus_inv_, cl_double);
//   set_kernel(d_ex_spikes_);
//   set_kernel(d_in_spikes_);
//   set_kernel_prim(time_index, cl_int);
  
//   return 0;
// }

// int
// nest::iaf_psc_alpha_cpu::set_static_deliver_kernel_args(cl::Kernel *kernel, int num_nodes, int batch_size)
// {
//   cl_int ret;
//   int arg_idx = 0;

//   set_kernel_prim(num_nodes, cl_int);
//   set_kernel_prim(event_size, cl_int);
//   set_kernel(d_connections_ptr);
//   set_kernel(d_connections);
//   set_kernel(d_connections_weight);
//   set_kernel_prim(batch_size, cl_int);
//   set_kernel(d_spike_src);
//   set_kernel(d_spike_multiplicity);
//   set_kernel(d_spike_pos);
//   set_kernel(d_ex_spikes_);
//   set_kernel(d_in_spikes_);
//   set_kernel_prim(time_index, cl_int);

//   return 0;
// }

// void
// nest::iaf_psc_alpha_cpu::execute_kernel(cl::Kernel *kernel, clContext_ *clCxt, size_t num_nodes)
// {
//   cl_int ret;
//   size_t local_work_size = 128;
//   size_t global_work_size = num_nodes % local_work_size == 0 ? num_nodes : (num_nodes / local_work_size + 1) * local_work_size;
  
//   ret = this->command_queue.enqueueNDRangeKernel(*(kernel), cl::NullRange, global_work_size, local_work_size);

//   if (ret != CL_SUCCESS)
//   {
//     printf("Failed to EnqueueNDRangeKernel. error code: %d\n", ret);
//     return;
//   }

//   //clFinish(this->command_queue);
// }

void
nest::iaf_psc_alpha_cpu::fill_spike_event_buffer( Event& e)
{
}

void
nest::iaf_psc_alpha_cpu::fill_event_buffer( SecondaryEvent& e)
{
  // GapJunctionEvent *gap_event = (GapJunctionEvent*)&e;
  // index sgid = gap_event->get_sender_gid() - 1;
  // double *h_sgid_event_buffer = h_event_buffer + event_size * sgid;

  // h_event_weight[sgid] = e.get_weight();

  // size_t i = 0;
  // std::vector< unsigned int >::iterator it = gap_event->begin();
  // // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  // while ( it != gap_event->end() )
  //   {
  //     h_sgid_event_buffer[ i ] =
  //     gap_event->get_weight() * gap_event->get_coeffvalue( it );
  //     ++i;
  //   }
}

// void
// nest::iaf_psc_alpha_cpu::fill_buffer_zero_double(clContext_ *clCxt, cl::Buffer &buffer, int size)
// {
//   cl_double zero = 0.0;
//   cl_int ret;

//   ret = this->command_queue.enqueueFillBuffer(buffer, zero, 0, size);
  
//   if (ret != CL_SUCCESS)
//     {
//       printf("Failed to EnqueueNDRangeKernel. error code: %d\n", ret);
//       return;
//     }
// }

// void
// nest::iaf_psc_alpha_cpu::fill_buffer_zero_uint(clContext_ *clCxt, cl::Buffer &buffer, int size)
// {
//   cl_uint zero = 0;
//   cl_int ret;

//   ret = this->command_queue.enqueueFillBuffer(buffer, zero, 0, size);
  
//   if (ret != CL_SUCCESS)
//     {
//       printf("Failed to enqueueFillBuffer. error code: %d\n", ret);
//       return;
//     }
// }


void
nest::iaf_psc_alpha_cpu::initialize()
{
  if (not is_ring_buffer_ready)
  {
    //if (initialize_command_queue())
    //  return;

    cout << "initialize: Setup Ring Buffers" << endl;
    event_size = kernel().connection_manager.get_min_delay()
      + kernel().connection_manager.get_max_delay();
    ring_buffer_size = (this->num_local_nodes - OFFSET) * event_size;
    cout << "RingBuffer Size: " << ring_buffer_size << " " << this->num_local_nodes << " X " << event_size << endl;
    h_currents_.resize(ring_buffer_size);// = new double[ring_buffer_size];
    h_ex_spikes_.resize(ring_buffer_size);// = new double[ring_buffer_size];
    h_in_spikes_.resize(ring_buffer_size);// = new double[ring_buffer_size];
    // h_currents__mark = new int[ring_buffer_size];
    // h_ex_spikes__mark = new int[ring_buffer_size];
    // h_in_spikes__mark = new int[ring_buffer_size];
    // h_currents__count = new int[total_num_nodes];
    // h_ex_spikes__count = new int[total_num_nodes];
    // h_in_spikes__count = new int[total_num_nodes];
    // h_currents__index = new int[ring_buffer_size];
    // h_ex_spikes__index = new int[ring_buffer_size];
    // h_in_spikes__index = new int[ring_buffer_size];

    // create(&gpu_context, &d_currents__buf, ring_buffer_size*sizeof(double));
    // create(&gpu_context, &d_ex_spikes__buf, ring_buffer_size*sizeof(double));
    // create(&gpu_context, &d_in_spikes__buf, ring_buffer_size*sizeof(double));

    //create(&gpu_context, &d_currents_, ring_buffer_size*sizeof(double));
    //create(&gpu_context, &d_ex_spikes_, ring_buffer_size*sizeof(double));
    //create(&gpu_context, &d_in_spikes_, ring_buffer_size*sizeof(double));

    // create(&gpu_context, &d_currents__count, total_num_nodes*sizeof(int));
    // create(&gpu_context, &d_ex_spikes__count, total_num_nodes*sizeof(int));
    // create(&gpu_context, &d_in_spikes__count, total_num_nodes*sizeof(int));
    // create(&gpu_context, &d_currents__index, ring_buffer_size*sizeof(int));
    // create(&gpu_context, &d_ex_spikes__index, ring_buffer_size*sizeof(int));
    // create(&gpu_context, &d_in_spikes__index, ring_buffer_size*sizeof(int));
    // fill_buffer_zero_double(&gpu_context, d_currents__buf, ring_buffer_size * sizeof(double));
    // fill_buffer_zero_double(&gpu_context, d_ex_spikes__buf, ring_buffer_size * sizeof(double));
    // fill_buffer_zero_double(&gpu_context, d_in_spikes__buf, ring_buffer_size * sizeof(double));

    //fill_buffer_zero_double(&gpu_context, d_currents_, ring_buffer_size * sizeof(double));
    //fill_buffer_zero_double(&gpu_context, d_ex_spikes_, ring_buffer_size * sizeof(double));
    //fill_buffer_zero_double(&gpu_context, d_in_spikes_, ring_buffer_size * sizeof(double));
    
    for(size_t i = 0; i < ring_buffer_size; i++) {
      h_currents_[i] = 0.0;
      h_ex_spikes_[i] = 0.0;
      h_in_spikes_[i] = 0.0;
    }

    is_ring_buffer_ready = true;
  }

  for (vector<vector<connection_info> >::iterator tgt_it = connections.begin();
       tgt_it != connections.end();
       tgt_it++)
    {
      graph_size += (*tgt_it).size();
    }

  cout << "Initialize Graph " << graph_size << endl;
  cout << "Connections " << connections.size() << " Num Nodes " << this->total_num_nodes << endl;

  h_connections_ptr.resize(this->total_num_nodes + 1); // = new int[this->total_num_nodes + 1];
  h_connections.resize(graph_size); // = new int[graph_size];
  h_connections_weight.resize(graph_size); // = new double[graph_size];
  //create(&gpu_context, &d_connections_ptr, (this->total_num_nodes + 1) * sizeof(int));
  //create(&gpu_context, &d_connections, graph_size * sizeof(int));
  //create(&gpu_context, &d_connections_weight, graph_size * sizeof(double));

  int i = 0;
  int node_id = 0;
  for (vector<vector<connection_info> >::iterator src_it = connections.begin();
       src_it != connections.end();
       src_it++, node_id++)
  {
      h_connections_ptr[node_id] = i;

      if ((*src_it).empty()) {
        cout << "Connection " << node_id << " empty" << endl;
      }

      for (vector<connection_info>::iterator tgt_it = (*src_it).begin();
       tgt_it != (*src_it).end();
       tgt_it++, i++)
    {
      connection_info& info = *tgt_it;
      h_connections[i] = info.tgt_id;
      h_connections_weight[i] = info.weight;
      // if (node_id == 68)
      //   cout << node_id << " " << info.tgt_id << " " << info.weight << endl;
    }
  }

  cout << "Last node_id " << node_id << endl;

  //for (size_t i = connections.size(); i <= this->total_num_nodes; i++)
  //  h_connections_ptr[i] = graph_size;
  h_connections_ptr[this->total_num_nodes] = graph_size;
  //synchronize();

  //upload(&gpu_context, h_connections_ptr, d_connections_ptr, (this->total_num_nodes + 1) * sizeof(int));
  //upload(&gpu_context, h_connections, d_connections, graph_size * sizeof(int));
  //upload(&gpu_context, h_connections_weight, d_connections_weight, graph_size * sizeof(double));

  //synchronize();

  //getKernel("update", "deliver_events_stdp_pl", "deliver_events", &gpu_context);

  //MOD
  //delete[] h_connections_ptr; h_connections_ptr = NULL;
  //delete[] h_connections; h_connections = NULL;
  //delete[] h_connections_weight; h_connections_weight = NULL;
}

void
nest::iaf_psc_alpha_cpu::clear_buffer()
{
  // for (unsigned int i = 0; i < ring_buffer_size; i++)
  // {
  //     h_currents_[i] = 0.0;
  //     h_ex_spikes_[i] = 0.0;
  //     h_in_spikes_[i] = 0.0;
  // }

// h_currents__mark[i] = 0;
// h_ex_spikes__mark[i] = 0;
// h_in_spikes__mark[i] = 0;
   //}
  //   h_event_buffer[i] = 0.0;
  
  //fill_buffer_zero_double(&gpu_context, d_ring_buffer, ring_buffer_size);
}

void
nest::iaf_psc_alpha_cpu::deliver_events()
{
  // update_type = 1 -- first update
  // update_type = 2 -- regular update (with multiplicity)

  size_t batch_size = list_spikes.size();
  //cout << "batch_size " << batch_size << endl;
  // getchar();

  //int thrd_id = kernel().vp_manager.get_thread_id();
  //#pragma omp critical
  //printf("[%d] Deliver Events - Batch Size: %d\n", thrd_id, batch_size);
  if (batch_size != 0)
  {
#ifdef PROFILING_SIZE
    #pragma omp critical
    printf("{P} Normal Batch: %d\n", batch_size);
#endif
    cout << "Normal Batch: " << batch_size << endl;

    int* h_spike_tgid = new int[batch_size];
    // h_t_spike = new double[batch_size];
    // h_dendritic_delay = new double[batch_size];
    double* h_weight_ = new double[batch_size];
    long* h_pos = new long[batch_size];
    double* h_t_lastspike = new double[batch_size];
    double* h_Kplus_ = new double[batch_size];
    int* h_conn_type_ = new int[batch_size];

    int* h_multiplicity;
    if (this->update_type == 2)
    {
      h_multiplicity = new int[batch_size];
    }

    //cout << "1" << endl;
    //create(&gpu_context, &d_spike_tgid, batch_size*sizeof(int));

    // create(&gpu_context, &d_t_spike, batch_size*sizeof(double));
    // create(&gpu_context, &d_dendritic_delay, batch_size*sizeof(double));
    
    //create(&gpu_context, &d_weight_, batch_size*sizeof(double));
    //create(&gpu_context, &d_pos, batch_size*sizeof(long));
    //create(&gpu_context, &d_t_lastspike, batch_size*sizeof(double));
    //create(&gpu_context, &d_Kplus_, batch_size*sizeof(double));
    //create(&gpu_context, &d_conn_type_, batch_size*sizeof(int));

    //if (this->update_type == 2)
    //{
    //  create(&gpu_context, &d_multiplicity, batch_size*sizeof(int));
    //}

    //synchronize();

    //cout << "2" << endl;
    int type_count = 0;
    int ind = 0;
    for (vector< synapse_info >::iterator it = list_spikes.begin(); it != list_spikes.end(); it++, ind++)
    {
      assert(ind < batch_size);

      synapse_info& entry = *it;
      h_spike_tgid[ind] = entry.target_node;
      // h_t_spike[ind] = entry.t_spike;
      // h_dendritic_delay[ind] = entry.dendritic_delay;
      h_weight_[ind] = entry.weight;
      
      h_conn_type_[ind] = entry.type;
      h_pos[ind] = entry.pos;
      if (this->update_type == 2)
      {
        h_multiplicity[ind] = entry.multiplicity;
      }

      if (entry.type == 2)
      {
        type_count++;
        h_Kplus_[ind] = entry.Kplus;
        h_t_lastspike[ind] = entry.last_t_spike;
      }
    }

    //cout << "3" << endl;
    //upload(&gpu_context, (void*)h_spike_tgid, d_spike_tgid, batch_size*sizeof(int));

    // upload(&gpu_context, (void*)h_t_spike, d_t_spike, batch_size*sizeof(double));
    // upload(&gpu_context, (void*)h_dendritic_delay, d_dendritic_delay, batch_size*sizeof(double));
    
    //upload(&gpu_context, (void*)h_weight_, d_weight_, batch_size*sizeof(double));
    //upload(&gpu_context, (void*)h_pos, d_pos, batch_size*sizeof(long));

    //if (this->update_type == 2)
    //{
    //  upload(&gpu_context, (void*)h_multiplicity, d_multiplicity, batch_size*sizeof(int));
    //}
  
    if (type_count > 0)
    {
      //upload(&gpu_context, (void*)h_t_lastspike, d_t_lastspike, batch_size*sizeof(double));
      //upload(&gpu_context, (void*)h_Kplus_, d_Kplus_, batch_size*sizeof(double));
      //upload(&gpu_context, (void*)h_conn_type_, d_conn_type_, batch_size*sizeof(int));
      conn_type = 0;
    }
    else
    {
      conn_type = 1;
    }
  
    //synchronize();

    //cout << "4 - type_count " << type_count << endl;
    //set_deliver_kernel_args(this->deliver_kernel, this->num_local_nodes, batch_size);
    //execute_kernel(this->deliver_kernel, &gpu_context, batch_size);

    for(size_t tid = 0; tid < batch_size; tid++) {
kernel_deliver_events_stdp_pl(
  tid,
  this->num_local_nodes - OFFSET,
  batch_size,
             event_size,
             // History
             h_history_ptr,
             h_history_Kminus_,
             h_history_t_,
             h_history_access_counter_,
             h_Kminus_,
             h_tau_minus_inv_,
             //
             h_spike_tgid,
             event_t_spike,
             event_dendritic_delay,
             h_weight_,
             h_Kplus_,
             h_conn_type_,
             h_t_lastspike,
             h_pos,
             h_multiplicity,
             event_multiplicity,
             this->update_type,
             conn_type,
             cp_lambda_,
             cp_mu_,
             cp_alpha_,
             cp_tau_plus_inv_,
             h_ex_spikes_,
             h_in_spikes_,
             time_index);
  }

    //synchronize();

    //cout << "5" << endl;
    cout << "conn_type " << conn_type << " type_count " << type_count << endl;
    int connection0 = 0, connection1 = 0, connection2 = 0;
    if (type_count > 0)
    {
      //download(&gpu_context, d_weight_, (void*)h_weight_, batch_size * sizeof(double));
      //download(&gpu_context, d_Kplus_, (void*)h_Kplus_, batch_size * sizeof(double));
  
      //synchronize();
      ind = 0;
      for (vector< synapse_info >::iterator it = list_spikes.begin(); it != list_spikes.end(); it++, ind++)
      {
        assert(ind < batch_size);

        synapse_info& entry = *it;
        if (entry.type == 1)
        {
          // STDPConnection< TargetIdentifierIndex > *connection = entry.connection;
          // connection->weight_ = h_weight_[ind];
          // connection->Kplus_ = h_Kplus_[ind];
          connection1 +=1;
        }
        else if (entry.type == 2)
        {
          STDPPLConnectionHom< TargetIdentifierIndex > *connection = entry.connection;
          connection->weight_ = h_weight_[ind];
          connection->Kplus_ = h_Kplus_[ind];

          //cout << "[" << ind << "] W" << h_weight_[ind] << endl;
          connection2 +=1;
        } else {
          connection0 +=1;
        }
      }
    }
    cout << "Updated Connection " << connection0 << " " << connection1 << " " << connection2 << endl;

    //cout << "6" << endl;
    delete[] h_spike_tgid;
    // delete[] h_t_spike;
    // delete[] h_dendritic_delay;
    delete[] h_weight_;
    delete[] h_pos;
    delete[] h_t_lastspike;
    if (this->update_type == 2)
    {
      delete[] h_multiplicity;
    }
    delete[] h_Kplus_;
    delete[] h_conn_type_;
  
    list_spikes.clear();
    //cout << "7" << endl;
  }
}

void
nest::iaf_psc_alpha_cpu::deliver_static_events()
{
  size_t static_batch_size = list_sgid.size();
  // cout << "static_batch_size " << static_batch_size << endl;
  // getchar();
  
  //int thrd_id = kernel().vp_manager.get_thread_id();
  //#pragma omp critical
  //printf("[%d] Static Batch: %d\n", thrd_id, static_batch_size);
  if (static_batch_size != 0)
  {
#ifdef PROFILING_SIZE
    #pragma omp critical
    printf("{P} Static Batch: %d\n", static_batch_size);
#endif
    cout << "Static Batch: " << static_batch_size << endl;

    int* h_spike_src = new int[static_batch_size];
    int* h_spike_multiplicity = new int[static_batch_size];
    long* h_spike_pos = new long[static_batch_size];

    //create(&gpu_context, &d_spike_src, static_batch_size*sizeof(int));
    //create(&gpu_context, &d_spike_multiplicity, static_batch_size*sizeof(int));
    //create(&gpu_context, &d_spike_pos, static_batch_size*sizeof(long));

    int i = 0;
    for (vector< SpikeEvent >::iterator it = list_sgid.begin(); it != list_sgid.end(); it++, i++)
    {
      SpikeEvent &e = *it;
      long pos = e.get_rel_delivery_steps(
                          kernel().simulation_manager.get_slice_origin() );
      h_spike_src[i] = e.get_sender_gid();
      h_spike_multiplicity[i] = e.get_multiplicity();
      h_spike_pos[i] = pos;
    }
    //upload(&gpu_context, (void*)h_spike_src, d_spike_src, static_batch_size*sizeof(int));
    //upload(&gpu_context, (void*)h_spike_multiplicity, d_spike_multiplicity, static_batch_size*sizeof(int));
    //upload(&gpu_context, (void*)h_spike_pos, d_spike_pos, static_batch_size*sizeof(long));

    //synchronize();
    //set_static_deliver_kernel_args(this->static_deliver_kernel, this->num_local_nodes, static_batch_size);
    //execute_kernel(this->static_deliver_kernel, &gpu_context, static_batch_size * 32 * 4);

    kernel_deliver_events(this->num_local_nodes - OFFSET,
           event_size,
           h_connections_ptr,
           h_connections,
           h_connections_weight,
           static_batch_size,
           h_spike_src,
           h_spike_multiplicity,
           h_spike_pos,
           h_ex_spikes_,
           h_in_spikes_,
           time_index);

    //synchronize();

    delete[] h_spike_src;
    delete[] h_spike_multiplicity;
    delete[] h_spike_pos;
    list_sgid.clear();
  }
}

void
nest::iaf_psc_alpha_cpu::copy_event_data(std::vector<Node *> nodes)
{
}

void
nest::iaf_psc_alpha_cpu::handle(index sgid, index tgid, double weight_)
{
  //assert(sgid < 11250);
  connection_info info;
  info.tgt_id = tgid;
  info.weight = weight_;
  connections[sgid].push_back(info);
}

void nest::iaf_psc_alpha_cpu::handle( SpikeEvent& e )
{};
void nest::iaf_psc_alpha_cpu::handle( CurrentEvent& e )
{};

void nest::iaf_psc_alpha_cpu::pre_deliver_event()
{
  this->pre_deliver_event(actualNodes);
}

void nest::iaf_psc_alpha_cpu::pre_deliver_event(const std::vector< Node* > &nodes)
{
  size_t len = nodes.size();
  
  // if (h_history_ptr == NULL)
  //   {
  //     h_history_ptr = new int[len + 1]{0};
  //     //create(&gpu_context, &d_history_ptr, (len + 1)*sizeof(int));
  //   }
  
  // if (h_Kminus_ == NULL)
  //   {
  //     h_Kminus_ = new double[len]{0};
  //     //create(&gpu_context, &d_Kminus_, len*sizeof(double));
  //   }
  
  // if (h_tau_minus_inv_ == NULL)
  //   {
  //     h_tau_minus_inv_ = new double[len]{0};
  //     //create(&gpu_context, &d_tau_minus_inv_, len*sizeof(double));
  //   }

  if (!isHistory) {
     h_history_ptr.resize(len + 1); // = new int[len + 1]{0};
     h_Kminus_.resize(len); //= new double[len]{0};
     h_tau_minus_inv_.resize(len); // = new double[len]{0};
     isHistory = true;
  }

  //synchronize();
  
  typedef std::deque< histentry > hist_queue;
  hist_queue nodes_history;
  
  int nodeid = 0;

  for (std::vector< Node* >::const_iterator it = nodes.begin(); it != nodes.end(); it++, nodeid++)
  {
    //if (nodeid >= this->num_local_nodes - OFFSET)
    //  break;

    nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*it;
    //nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)nodes.at(nodeid);
    
    h_history_ptr[nodeid] = nodes_history.size();
    //history_size = nodes_history.size();
    hist_queue &h_ = node->get_all_history();
    hist_queue::iterator last_it = nodes_history.end();
    nodes_history.insert(last_it, h_.begin(), h_.end());
    h_Kminus_[nodeid] = node->get_Kminus_();
    h_tau_minus_inv_[nodeid] = node->get_tau_minus_inv_();
  }
  
  history_size = nodes_history.size();

  cout << "pre_deliver_event: Len " << len << " History " << history_size << endl;
  
  //for (size_t nodeid = this->num_local_nodes; nodeid <= len; nodeid++)
  //  h_history_ptr[nodeid] = history_size;
  h_history_ptr[len] = history_size;

  //upload(&gpu_context, (void*)h_history_ptr, d_history_ptr, (len + 1) * sizeof(int));
  //upload(&gpu_context, (void*)h_Kminus_, d_Kminus_, len*sizeof(double));
  //upload(&gpu_context, (void*)h_tau_minus_inv_, d_tau_minus_inv_, len*sizeof(double));

  //synchronize();
  
  h_history_Kminus_.resize(history_size); // = new double[history_size];
  h_history_t_.resize(history_size); // = new double[history_size];
  h_history_access_counter_.resize(history_size); // = new int[history_size];

  //create(&gpu_context, &d_history_Kminus_, history_size*sizeof(double));
  //create(&gpu_context, &d_history_t_, history_size*sizeof(double));
  //create(&gpu_context, &d_history_access_counter_, history_size*sizeof(int));

  //synchronize();
  
  int index = 0;
  for (hist_queue::iterator hist_it = nodes_history.begin();
       hist_it != nodes_history.end();
       hist_it++, index++)
    {
      histentry& entry = *hist_it;
      h_history_Kminus_[index] = entry.Kminus_;
      h_history_t_[index] = entry.t_;
      h_history_access_counter_[index] = entry.access_counter_;
    }

  //upload(&gpu_context, (void*)h_history_Kminus_, d_history_Kminus_, history_size*sizeof(double));
  //upload(&gpu_context, (void*)h_history_t_, d_history_t_, history_size*sizeof(double));
  //upload(&gpu_context, (void*)h_history_access_counter_, d_history_access_counter_, history_size*sizeof(int));
  //fill_buffer_zero_uint(&gpu_context, d_history_access_counter_, history_size*sizeof(int));      

  //synchronize();
}

void nest::iaf_psc_alpha_cpu::post_deliver_event()
{
  this->post_deliver_event(actualNodes);
}

void nest::iaf_psc_alpha_cpu::post_deliver_event(const std::vector< Node* > &nodes)
{
  typedef std::deque< histentry > hist_queue;
  //download(&gpu_context, d_history_access_counter_, (void*)h_history_access_counter_, history_size * sizeof(int));

  //synchronize();
  
  int hist_it = 0;
  int nodeid = 0;
  for (std::vector< Node* >::const_iterator it = nodes.begin(); it != nodes.end();
       it++, nodeid++)
  {
    //if (nodeid >= this->num_local_nodes - OFFSET)
    //  break;

    nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*it;
    hist_queue h_ = node->get_all_history();
    for (hist_queue::iterator h_it = h_.begin(); h_it != h_.end(); h_it++, hist_it++)
    {
      h_it->access_counter_ = h_history_access_counter_[hist_it];
    }
  }

  //delete[] h_history_Kminus_; h_history_Kminus_ = NULL;
  //delete[] h_history_t_; h_history_t_ = NULL;
  //delete[] h_history_access_counter_; h_history_access_counter_ = NULL;
}

void nest::iaf_psc_alpha_cpu::handle(Event& ev, double last_t_spike, const CommonSynapseProperties *csp, void *conn, int conn_type)
{
  SpikeEvent *e = (SpikeEvent *)&ev;

  // STDPPLConnectionHom::send
  if (conn_type == 2)
    {
      const STDPPLHomCommonProperties *cp = (STDPPLHomCommonProperties *)csp;
      int sgid = e->get_sender_gid();
      int tgid = e->get_receiver().get_thread_lid();

      cp_lambda_ = cp->lambda_;
      cp_mu_ = cp->mu_;
      cp_alpha_ = cp->alpha_;
      cp_tau_plus_inv_ = cp->tau_plus_inv_;

      double t_spike = e->get_stamp().get_ms();

      STDPPLConnectionHom< TargetIdentifierIndex > *connection = (STDPPLConnectionHom< TargetIdentifierIndex > *) conn;

      synapse_info conn_info;
      conn_info.source_node = sgid;
      conn_info.target_node = tgid;
      conn_info.connection = connection;
      //conn_info.t_spike = t_spike;
      conn_info.last_t_spike = last_t_spike;
      //conn_info.dendritic_delay = connection->get_delay();
      conn_info.weight = connection->weight_;
      conn_info.Kplus = connection->Kplus_;
      conn_info.pos = e->get_rel_delivery_steps(
                        kernel().simulation_manager.get_slice_origin() );
      conn_info.multiplicity = e->get_multiplicity();
      conn_info.type = conn_type;
      list_spikes.push_back(conn_info);

      event_t_spike = t_spike;
      event_dendritic_delay = connection->get_delay();

      event_multiplicity = e->get_multiplicity();
    }
  else if (conn_type == 1)
    {
      int sgid = e->get_sender_gid();
      int tgid = e->get_receiver().get_thread_lid();

      synapse_info conn_info;
      conn_info.source_node = sgid;
      conn_info.target_node = tgid;
      conn_info.weight = e->get_weight();
      conn_info.pos = e->get_rel_delivery_steps(
                              kernel().simulation_manager.get_slice_origin() );
      conn_info.multiplicity = e->get_multiplicity();
      conn_info.type = conn_type;
      list_spikes.push_back(conn_info);

      event_multiplicity = e->get_multiplicity();
    }
}

void nest::iaf_psc_alpha_cpu::insert_static_event(SpikeEvent& e)
{
  //SpikeEvent *e = (SpikeEvent *)&ev;

  /*int sgid = e.get_sender_gid();
  int tgid = e.get_receiver().get_thread_lid();

  synapse_info conn_info;
  conn_info.source_node = sgid;
  conn_info.target_node = tgid;
  conn_info.pos = e.get_rel_delivery_steps(
                       kernel().simulation_manager.get_slice_origin() );
  conn_info.weight = e.get_weight();
  conn_info.multiplicity = e.get_multiplicity();
  conn_info.type = 1;
  list_spikes.push_back(conn_info);*/
  list_sgid.push_back(e);
}

// iaf_psc_alpha::handle(SpikeEvent& e)
void nest::iaf_psc_alpha_cpu::insert_event(SpikeEvent& e)
{
  int sgid = e.get_sender_gid();
  int tgid = e.get_receiver().get_thread_lid();

  synapse_info conn_info;
  conn_info.source_node = sgid;
  conn_info.target_node = tgid;
  conn_info.pos = e.get_rel_delivery_steps(
                       kernel().simulation_manager.get_slice_origin() );
  conn_info.weight = e.get_weight();
  conn_info.multiplicity = e.get_multiplicity();
  conn_info.type = 1;
  list_spikes.push_back(conn_info);
}

void
nest::iaf_psc_alpha_cpu::advance_time()
{
  time_index = (time_index + kernel().connection_manager.get_min_delay()) % event_size;
}

inline
bool nest::iaf_psc_alpha_cpu::isLocalNode(nest::index i) const {
  return i < this->num_local_nodes - OFFSET;
}

#include <fstream>
void writeNodes(nest::thread t, int count, const std::vector< nest::Node* > &nodes)
{
  std::ofstream myfile;
  myfile.open ("C" + std::to_string(t) + "-" + std::to_string(count) + ".txt");

  for (std::vector<nest::Node *>::const_iterator nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
  {
    nest::iaf_psc_alpha* node = (nest::iaf_psc_alpha*)*nodeIt;
    myfile << node->S_.y3_ << endl;
  }

  myfile.close();
}

void writeNodes(nest::thread t, int count, double* nodes, size_t num_nodes)
{
  std::ofstream myfile;
  myfile.open ("G" + std::to_string(t) + "-" + std::to_string(count) + ".txt");

  for(size_t i = 0; i < num_nodes; i++) 
  {
    myfile << nodes[i] << endl;
  }

  myfile.close();
}
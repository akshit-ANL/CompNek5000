#include <catalyst.h>
#include <conduit.h>

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace
{
bool initialized = false;
int rank_id = 0;

void check_status(int status, const char* operation)
{
  if (status != catalyst_status_ok)
  {
    std::cerr << "Catalyst 2: " << operation << " failed with status "
              << status << std::endl;
    std::abort();
  }
}
}

extern "C" void nek_catalyst_initialize_(const int* rank)
{
  rank_id = *rank;
  const char* script = std::getenv("NEK_CATALYST_SCRIPT");
  if (!script || !*script)
    script = "catalyst_pipeline.py";

  conduit_node* params = conduit_node_create();
  conduit_node_set_path_char8_str(
    params, "catalyst/scripts/nek/filename", script);
  check_status(catalyst_initialize(params), "initialize");
  conduit_node_destroy(params);
  initialized = true;

  if (rank_id == 0)
    std::cout << "Catalyst 2 initialized with pipeline " << script << std::endl;
}

extern "C" void nek_catalyst_execute_(
  const double* x, const double* y, const double* z,
  const double* vx, const double* vy, const double* vz,
  const double* temperature, const int* lx1, const int* ly1,
  const int* lz1, const int* nelt, const int* ndim,
  const int* has_temperature, const int* step, const double* time)
{
  if (!initialized)
    return;

  const conduit_index_t nx = *lx1;
  const conduit_index_t ny = *ly1;
  const conduit_index_t nz = (*ndim == 3) ? *lz1 : 1;
  const conduit_index_t points_per_element = nx * ny * nz;
  const conduit_index_t number_of_points = points_per_element * *nelt;
  const conduit_index_t cells_per_element =
    (nx - 1) * (ny - 1) * ((*ndim == 3) ? (nz - 1) : 1);
  const conduit_index_t vertices_per_cell = (*ndim == 3) ? 8 : 4;

  std::vector<std::int64_t> connectivity;
  connectivity.reserve(static_cast<std::size_t>(
    cells_per_element * *nelt * vertices_per_cell));
  for (int element = 0; element < *nelt; ++element)
  {
    const std::int64_t offset =
      static_cast<std::int64_t>(element) * points_per_element;
    for (conduit_index_t k = 0; k < ((*ndim == 3) ? nz - 1 : 1); ++k)
      for (conduit_index_t j = 0; j < ny - 1; ++j)
        for (conduit_index_t i = 0; i < nx - 1; ++i)
        {
          const std::int64_t p0 = offset + i + nx * (j + ny * k);
          const std::int64_t p1 = p0 + 1;
          const std::int64_t p3 = p0 + nx;
          const std::int64_t p2 = p3 + 1;
          connectivity.insert(connectivity.end(), {p0, p1, p2, p3});
          if (*ndim == 3)
          {
            const std::int64_t plane = nx * ny;
            connectivity.insert(connectivity.end(),
              {p0 + plane, p1 + plane, p2 + plane, p3 + plane});
          }
        }
  }

  conduit_node* params = conduit_node_create();
  conduit_node_set_path_int64(params, "catalyst/state/timestep", *step);
  conduit_node_set_path_float64(params, "catalyst/state/time", *time);
  conduit_node_set_path_char8_str(params, "catalyst/channels/mesh/type", "mesh");

  const std::string mesh = "catalyst/channels/mesh/data/";
  conduit_node_set_path_char8_str(params,
    (mesh + "coordsets/coords/type").c_str(), "explicit");
  conduit_node_set_path_external_float64_ptr(params,
    (mesh + "coordsets/coords/values/x").c_str(), const_cast<double*>(x),
    number_of_points);
  conduit_node_set_path_external_float64_ptr(params,
    (mesh + "coordsets/coords/values/y").c_str(), const_cast<double*>(y),
    number_of_points);
  conduit_node_set_path_external_float64_ptr(params,
    (mesh + "coordsets/coords/values/z").c_str(), const_cast<double*>(z),
    number_of_points);
  conduit_node_set_path_char8_str(params,
    (mesh + "topologies/mesh/type").c_str(), "unstructured");
  conduit_node_set_path_char8_str(params,
    (mesh + "topologies/mesh/coordset").c_str(), "coords");
  conduit_node_set_path_char8_str(params,
    (mesh + "topologies/mesh/elements/shape").c_str(),
    (*ndim == 3) ? "hex" : "quad");
  conduit_node_set_path_external_int64_ptr(params,
    (mesh + "topologies/mesh/elements/connectivity").c_str(),
    connectivity.data(), static_cast<conduit_index_t>(connectivity.size()));

  conduit_node_set_path_char8_str(params,
    (mesh + "fields/velocity/association").c_str(), "vertex");
  conduit_node_set_path_char8_str(params,
    (mesh + "fields/velocity/topology").c_str(), "mesh");
  conduit_node_set_path_external_float64_ptr(params,
    (mesh + "fields/velocity/values/u").c_str(), const_cast<double*>(vx),
    number_of_points);
  conduit_node_set_path_external_float64_ptr(params,
    (mesh + "fields/velocity/values/v").c_str(), const_cast<double*>(vy),
    number_of_points);
  conduit_node_set_path_external_float64_ptr(params,
    (mesh + "fields/velocity/values/w").c_str(), const_cast<double*>(vz),
    number_of_points);

  if (*has_temperature)
  {
    conduit_node_set_path_char8_str(params,
      (mesh + "fields/temperature/association").c_str(), "vertex");
    conduit_node_set_path_char8_str(params,
      (mesh + "fields/temperature/topology").c_str(), "mesh");
    conduit_node_set_path_external_float64_ptr(params,
      (mesh + "fields/temperature/values").c_str(),
      const_cast<double*>(temperature), number_of_points);
  }

  check_status(catalyst_execute(params), "execute");
  conduit_node_destroy(params);
}

extern "C" void nek_catalyst_finalize_()
{
  if (!initialized)
    return;
  conduit_node* params = conduit_node_create();
  check_status(catalyst_finalize(params), "finalize");
  conduit_node_destroy(params);
  initialized = false;
  if (rank_id == 0)
    std::cout << "Catalyst 2 finalized" << std::endl;
}

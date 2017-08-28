#ifndef HELPER_H
#define HELPER_H

#include "mfem.hpp"

#include <memory>
#include <string>
#include <sstream>

struct Simulation {
  std::shared_ptr<mfem::Mesh> mesh;
  mfem::GridFunction sol;
  mfem::GridFunction sol2;
};

inline Simulation load_sim(const std::string &str, const std::string separator) {
  Simulation sim;
  auto sep_pos = str.find(separator);

  {
    std::istringstream meshreader(std::string(str, 0, sep_pos));
    sim.mesh.reset(new mfem::Mesh(meshreader));
  }

  auto sep2_pos = str.find(separator, sep_pos + separator.size());
  {
    std::istringstream gridreader(
        std::string(str, sep_pos + separator.size(), sep2_pos));
    sim.sol = mfem::GridFunction(sim.mesh.get(), gridreader);
  }

  if (sep2_pos != std::string::npos) {
    std::istringstream gridreader(
        std::string(str, sep2_pos + separator.size(), std::string::npos));
    sim.sol2 = mfem::GridFunction(sim.mesh.get(), gridreader);
  }

  return sim;
}

#endif // HELPER_H

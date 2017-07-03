#ifndef HELPER_H
#define HELPER_H


struct Simulation {
  std::shared_ptr<mfem::Mesh> mesh;
  mfem::GridFunction sol;
};

inline Simulation load_sim(const std::string &str, const std::string separator) {
  Simulation sim;
  auto sep_pos = str.find(separator);

  std::istringstream meshreader(std::string(str, 0, sep_pos));
  sim.mesh.reset(new mfem::Mesh(meshreader));

  std::istringstream gridreader(
      std::string(str, sep_pos + separator.size(), std::string::npos));
  sim.sol = mfem::GridFunction(sim.mesh.get(), gridreader);

  return sim;
}

#endif // HELPER_H

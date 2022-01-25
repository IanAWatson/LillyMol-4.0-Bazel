// 3d substructure searching on SMU

#include "riegeli/records/record_reader.h"
#include <tensorflow/core/platform/default/posix_file_system.h>
#include "tensorflow/core/platform/env.h"

#include "Molecule_Lib/substructure.h"

namespace smu_3d_substructure_search {
int
substructure_search(std::string& fname,
                    Substructure_Query& query,
                    std::ostream& output) {
  Env* env = Env::Default();
//tensorflow::io::RecordReader reader(fname, 
  std::unique_ptr<RandomAccessFile> read_file;
  // Read it back with the RecordReader.
  env->NewRandomAccessFile(fname, &read_file);
}

int main(int argc, char** argv) {
}
}  // namespace smu_3d_substructure_search

int
main(int argc, char ** argv) {
  smu_3d_substructure_search::main();
  return 0;
}

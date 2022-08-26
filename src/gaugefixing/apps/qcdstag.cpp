#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/access_pattern/GpuPattern.hxx"
#include "../../lattice/access_pattern/StandardPattern.hxx"
#include "../../lattice/datatype/datatypes.h"
#include "include/c-lime/lime.h"
#include "include/c-lime/lime_config.h"
#include "include/c-lime/lime_fixed_types.h"
#include "string.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

bool readQCDSTAG(SiteCoord<4, FULL_SPLIT> s, const char *file_name,
                 const short SIZE[4], Real *U) {

  int data_size = 4 * SIZE[0] * SIZE[1] * SIZE[2] * SIZE[3] * 18;

  std::vector<double> v(data_size);

  std::ifstream stream(file_name);

  if (!stream.read((char *)&v[0], data_size * sizeof(double))) {
    std::cout << "readQCDSTAG error: " << file_name << std::endl;
    stream.close();
    return false;
  }

  //SiteCoord<4, NO_SPLIT> s1(SIZE);

  int unique;
  int mu1;
  for (int t = 0; t < SIZE[0]; t++) {
    for (int z = 0; z < SIZE[3]; z++) {
      for (int y = 0; y < SIZE[2]; y++) {
        for (int x = 0; x < SIZE[1]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            if (mu == 3)
              mu1 = 0;
            else
              mu1 = mu + 1;
            for (int j = 0; j < 3; j++) {
              for (int k = 0; k < 3; k++) {
                for (int c = 0; c < 2; c++) {
                  /*unique = (t * SIZE[1] * SIZE[2] * SIZE[3] +
                            x * SIZE[3] * SIZE[2] + y * SIZE[3] + z) *
                               4 * 18 +
                           mu1 * 18 + j * 6 + k * 2 + c;*/
                  s.site[0] = t;
                  s.site[1] = x;
                  s.site[2] = y;
                  s.site[3] = z;
                  //unique = StandardPattern<SiteCoord<4, NO_SPLIT>, 4,
                                           //3>::getUniqueIndex(s1, mu1, j, k, c);
                  /*U[GpuPattern<SiteCoord<4, FULL_SPLIT>, 4,
                               3>::getIndexByUnique(unique, s.size)] =
                      v[(mu1 * SIZE[0] * SIZE[1] * SIZE[2] * SIZE[3] +
                         t * SIZE[1] * SIZE[2] * SIZE[3] +
                         z * SIZE[1] * SIZE[2] + y * SIZE[1] + x) *
                            18 +
                        j * 6 + k * 2 + c];*/
                  U[GpuPattern<SiteCoord<4, FULL_SPLIT>, 4, 3>::getIndex(s, mu, j, k, (bool)c)] =
                      v[(mu * SIZE[0] * SIZE[1] * SIZE[2] * SIZE[3] +
                         t * SIZE[1] * SIZE[2] * SIZE[3] +
                         z * SIZE[1] * SIZE[2] + y * SIZE[1] + x) * 18 +
                         j * 6 + k * 2 + c];
                  /*unique = (t * SIZE[1] * SIZE[2] * SIZE[3] +
                            x * SIZE[3] * SIZE[2] + y * SIZE[3] + z) *
                               4 * 18 +
                           mu1 * 18 + j * 6 + k * 2 + c;
                  U[GpuPattern<SiteCoord<4, FULL_SPLIT>, 4,
                               3>::getIndexByUnique(unique, s.size)] =
                      v[(mu1 * SIZE[0] * SIZE[1] * SIZE[2] * SIZE[3] +
                         t * SIZE[1] * SIZE[2] * SIZE[3] +
                         z * SIZE[1] * SIZE[2] + y * SIZE[1] + x) *
                            18 +
                        j * 6 + k * 2 + c];*/
                }
              }
            }
          }
        }
      }
    }
  }

  stream.close();

  return true;
}

bool writeQCDSTAG(SiteCoord<4, FULL_SPLIT> s, const char *output_name,
                  const short SIZE[4], Real *U) {

  int data_size = 4 * SIZE[0] * SIZE[1] * SIZE[2] * SIZE[3] * 18;

  std::vector<double> v(data_size);

  std::ofstream stream(output_name);

  int unique;
  int place;
  int mu1;
  for (int t = 0; t < SIZE[0]; t++) {
    for (int z = 0; z < SIZE[3]; z++) {
      for (int y = 0; y < SIZE[2]; y++) {
        for (int x = 0; x < SIZE[1]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            if (mu == 3)
              mu1 = 0;
            else
              mu1 = mu + 1;
            for (int j = 0; j < 3; j++) {
              for (int k = 0; k < 3; k++) {
                for (int c = 0; c < 2; c++) {
                  unique = (t * SIZE[1] * SIZE[2] * SIZE[3] +
                            x * SIZE[3] * SIZE[2] + y * SIZE[3] + z) *
                               4 * 18 +
                           mu1 * 18 + j * 6 + k * 2 + c;
                  place = (mu1 * SIZE[0] * SIZE[1] * SIZE[2] * SIZE[3] +
                           t * SIZE[1] * SIZE[2] * SIZE[3] +
                           z * SIZE[1] * SIZE[2] + y * SIZE[1] + x) *
                              18 +
                          j * 6 + k * 2 + c;
                  v[place] = U[GpuPattern<SiteCoord<4, FULL_SPLIT>, 4,
                                          3>::getIndexByUnique(unique, s.size)];
                }
              }
            }
          }
        }
      }
    }
  }

  if (!stream.write((char *)&v[0], (data_size) * sizeof(double))) {
    std::cout << "writeQCDSTAG error: " << output_name << std::endl;
    stream.close();
    return false;
  } else {
    stream.close();
    return true;
  }
}

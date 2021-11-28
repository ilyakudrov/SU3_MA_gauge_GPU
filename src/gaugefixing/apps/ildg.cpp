#include "include/c-lime/lime.h"
#include "include/c-lime/lime_config.h"
#include "include/c-lime/lime_fixed_types.h"
#include <cstdlib>
#include "string.h"
#include "../../lattice/datatype/datatypes.h"
#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/access_pattern/GpuPattern.hxx"
#include <iostream>
#include <string>
#include <sstream>

double reverseValue(const char *data) {
  double result;

  char *dest = (char *)&result;

  for (int i = 0; i < sizeof(double); i++) {
    dest[i] = data[sizeof(double) - i - 1];
  }
  return result;
}

void readILDG(SiteCoord<4,FULL_SPLIT> s, const char *file_name, const short SIZE[4], Real *U) {

  FILE *fp;
  fp = fopen(file_name, "r");

  LimeReader *reader;
  reader = limeCreateReader(fp);

  int status;
  char *lime_type;
  n_uint64_t nbytes;
  double *v;
  while ((status = limeReaderNextRecord(reader)) != LIME_EOF) {

    if (status != LIME_SUCCESS) {
      fprintf(stderr, "limeReaderNextRecord returned status = %d\n", status);
    }

    lime_type = limeReaderType(reader);
    nbytes = limeReaderBytes(reader);

    if (strcmp(lime_type, "ildg-binary-data") == 0) {

      v = (double *)malloc(nbytes);

      status = limeReaderReadData((char *)&v[0], &nbytes, reader);

      for (int i = 0; i < nbytes / sizeof(double); i++) {
        v[i] = reverseValue((char *)&v[i]);
      }

      if (status != LIME_SUCCESS) {
        fprintf(stderr, "limeReaderReadData returned status = %d\n", status);
      }

      int unique;
      int mu1;
      for (int t = 0; t < SIZE[0]; t++) {
        for (int z = 0; z < SIZE[3]; z++) {
          for (int y = 0; y < SIZE[2]; y++) {
            for (int x = 0; x < SIZE[1]; x++) {
              for (int mu = 0; mu < 4; mu++) {
                if(mu == 3) mu1 = 0;
                else mu1 = mu + 1;
                for (int j = 0; j < 3; j++) {
                  for (int k = 0; k < 3; k++) {
                    for(int c = 0;c < 2;c++){
                      unique = (t * SIZE[1] *SIZE[2] *SIZE[3] + x * SIZE[3] *SIZE[2] + y * SIZE[3] + z) * 4 * 18 + mu1 * 18 + j * 6 + k * 2 + c;
                      U[GpuPattern< SiteCoord<4,FULL_SPLIT>,4,3>::getIndexByUnique( unique , s.size )] =
                       v[(t * SIZE[1] *SIZE[2] *SIZE[3] + z * SIZE[1] *SIZE[2] + y * SIZE[1] + x) * 4 * 18 + mu * 18 + j * 6 + k * 2 + c];
                    }
                  }
                }
              }
            }
          }
        }
      }
      free(v);
    }
  }
  fclose(fp);
}

void writeILDG(SiteCoord<4,FULL_SPLIT> s, const char *file_name, const char *output_name, const short SIZE[4], Real *U, int steps) {

  // reader
  FILE *fp;
  fp = fopen(file_name, "r");

  LimeReader *reader;
  reader = limeCreateReader(fp);

  int status;
  char *lime_type;
  int bytes_pad;
  int MB_flag, ME_flag;
  n_uint64_t nbytes;
  double *v;
  char *data;

  // writer
  FILE *fp_out;
  fp_out = fopen(output_name, "w");

  LimeWriter *writer;
  writer = limeCreateWriter(fp_out);

  LimeRecordHeader *h;

  while ((status = limeReaderNextRecord(reader)) != LIME_EOF) {

    if (status != LIME_SUCCESS) {
      fprintf(stderr, "limeReaderNextRecord returned status = %d\n", status);
    }

    lime_type = limeReaderType(reader);
    nbytes = limeReaderBytes(reader);
    bytes_pad = limeReaderPadBytes(reader);
    MB_flag = limeReaderMBFlag(reader);
    ME_flag = limeReaderMEFlag(reader);

    // printf("\n\n");
    // printf("Type:           %s\n", lime_type);
    // printf("Data Length:    %ld\n", nbytes);
    // printf("Padding Length: %d\n", bytes_pad);
    // printf("MB flag:        %d\n", MB_flag);
    // printf("ME flag:        %d\n", ME_flag);

    if (strcmp(lime_type, "ildg-binary-data") == 0) {

      v = (double *)malloc(nbytes);

      // for(int i = 0; i < nbytes/sizeof(double);i++){
      //   v[i] = (double)U[i];
      //   v[i] = reverseValue((char *)&v[i]);
      // }

      int unique;
      int place;
      int mu1;
      for (int t = 0; t < SIZE[0]; t++) {
        for (int z = 0; z < SIZE[3]; z++) {
          for (int y = 0; y < SIZE[2]; y++) {
            for (int x = 0; x < SIZE[1]; x++) {
              for (int mu = 0; mu < 4; mu++) {
                if(mu == 3) mu1 = 0;
                else mu1 = mu + 1;
                for (int j = 0; j < 3; j++) {
                  for (int k = 0; k < 3; k++) {
                    for(int c = 0;c < 2;c++){
                      unique = (t * SIZE[1] *SIZE[2] *SIZE[3] + x * SIZE[3] *SIZE[2] + y * SIZE[3] + z) * 4 * 18 + mu1 * 18 + j * 6 + k * 2 + c;
                      place = (t * SIZE[1] *SIZE[2] *SIZE[3] + z * SIZE[1] *SIZE[2] + y * SIZE[1] + x) * 4 * 18 + mu * 18 + j * 6 + k * 2 + c;
                      v[place] = U[GpuPattern< SiteCoord<4,FULL_SPLIT>,4,3>::getIndexByUnique( unique , s.size )];
                      v[place] = reverseValue((char *)&v[place]);
                    }
                  }
                }
              }
            }
          }
        }
      }

      h = limeCreateHeader(MB_flag, ME_flag, lime_type, nbytes);

      status = limeWriteRecordHeader(h, writer);

      if (status < 0) {
        fprintf(stderr, "LIME write header error %d\n", status);
      }

      status = limeWriteRecordData(v, &nbytes, writer);
      if (status != LIME_SUCCESS)
        fprintf(stderr, "LIME write error %d\n", status);

      free(v);
    } else {

      data = (char *)malloc(nbytes);
      status = limeReaderReadData((char *)data, &nbytes, reader);
      if (status != LIME_SUCCESS)
        fprintf(stderr, "LIME read error %d\n", status);

      // std::cout << "data: " << data << std::endl;

      std::string data_new;

      if (strcmp(lime_type, "xlf-info") == 0){
	std::stringstream filename(std::stringstream::out);
	filename << std::string(data) + " SA steps "<<steps;
        data_new = filename.str();
      }
      else
        data_new = std::string(data);

      // std::cout << "data new: " << data_new << std::endl;

      nbytes = data_new.length();

      h = limeCreateHeader(MB_flag, ME_flag, lime_type, nbytes);

      status = limeWriteRecordHeader(h, writer);

      if (status < 0) {
        fprintf(stderr, "LIME write header error %d\n", status);
      }

      status = limeWriteRecordData((void*)data_new.c_str(), &nbytes, writer);
      if (status != LIME_SUCCESS)
        fprintf(stderr, "LIME write error %d\n", status);

      free(data);
    }
  }

  fclose(fp);
  fclose(fp_out);
}

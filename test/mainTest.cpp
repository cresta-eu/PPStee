/*! \file mainTest.cpp
 * \brief Main test file.
 *
 * \date 14.05.2013
 * \author matu_gr
 *
 * \copyright Copyright 2013 German Aerospace Center (http://www.DLR.de)\n\n
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <algorithm>
#include <vector>

#include <mpi.h>

#include "gtest/gtest.h"

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running main() from mainTest.cpp" << std::endl;

  testing::InitGoogleTest(&argc, argv);
  
  // MPI init
  MPI_Init(&argc, &argv);
  int mpiMe, mpiN;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiMe);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiN);
  int allowedThreads[] = {2, 3, 4, 5, 6, 12, 113};
  std::vector<int> vecAllowedThreads(allowedThreads, allowedThreads + sizeof(allowedThreads) / sizeof(int));
  if (std::find(vecAllowedThreads.begin(), vecAllowedThreads.end(), mpiN) == vecAllowedThreads.end()) {
    std::cout << "Start with ";
    for (std::vector<int>::iterator it=vecAllowedThreads.begin(); it!=vecAllowedThreads.end(); ++it) { std::cout << (it==vecAllowedThreads.begin()?"":", ") << *it; }
    std::cout << "MPI threads! Abort." << std::endl;
    return -1;
  }

  int retTest = RUN_ALL_TESTS();

  // MPI finalize
  MPI_Finalize();

  // exit
  return retTest;
}

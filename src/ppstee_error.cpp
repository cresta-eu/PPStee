/*! \file ppstee_error.cpp
 * \brief Error handling for PPStee.
 *
 * \date 30.01.2013
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

#include "ppstee_error.hpp"

string PPSteeError::getErrorInformation(PPSTEE_ERROR err) {
    switch (err) {
    case PPSTEE_ERROR_SUCCESS: return "Everything's fine.";
    case PPSTEE_ERROR_GRAPHDATAINCONSISTENT: return "Graph data is inconsistent.";
    case PPSTEE_ERROR_MATCHINGERROR: return "Graph and weights do not match.";
    case PPSTEE_ERROR_GRAPHUNDETERMINED: return "Unable to perform partitioning: graph not submitted or undetermined.";
    case PPSTEE_ERROR_PARTITIONERNOTINSTALLED: return "Unable to perform partitioning: partitioner is not installed.";
    case PPSTEE_ERROR_PARMETIS: return "ParMETIS call erroneous.";
    case PPSTEE_ERROR_PTSCOTCH: return "PTScotch call erroneous.";
    case PPSTEE_ERROR_UNKNOWNERROR: return "Something went wrong, but any further information is not available.";
    }
    return "Specified error code unknown!";
  }

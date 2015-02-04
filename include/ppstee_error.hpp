/*! \file ppstee_error.hpp
 * \brief Error handling for PPStee.
 *
 * \date 18.09.2012
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

#ifndef PPSTEE_ERROR_HPP
#define PPSTEE_ERROR_HPP

#include<string>


//! PPStee error enumeration.
/*!
 * Most PPStee calls provide error codes enumerated here. Use getErrorInformation() for details.
 */
typedef enum ENUM_PPSTEE_ERRORTYPE {
  PPSTEE_ERROR_SUCCESS,    //!< PPSTEE_ERROR_SUCCESS
  PPSTEE_ERROR_GRAPHDATAINCONSISTENT,   //!< PPSTEE_ERROR_GRAPHDATAINCONSISTENT
  PPSTEE_ERROR_MATCHINGERROR,  //!< PPSTEE_ERROR_MATCHINGERROR
  PPSTEE_ERROR_GRAPHUNDETERMINED, //!< PPSTEE_ERROR_GRAPHUNDETERMINED
  PPSTEE_ERROR_PARTITIONERNOTINSTALLED, //!< PPSTEE_ERROR_PARTITIONERNOTINSTALLED
  PPSTEE_ERROR_PARMETIS, //!< PPSTEE_ERROR_PARMETIS
  PPSTEE_ERROR_PTSCOTCH, //!< PPSTEE_ERROR_PTSCOTCH
  PPSTEE_ERROR_ZOLTAN, //!< PPSTEE_ERROR_ZOLTAN
  PPSTEE_ERROR_UNKNOWNERROR//!< PPSTEE_ERROR_UNKNOWNERROR
} PPSTEE_ERROR;


#include "ppstee.hpp"
#include "ppstee_objects.hpp"


using std::string;

//! PPStee error class
/*!
 * Essentially, this class provides error information for error codes.
 *
 * Usage:
 * string err_info = PPSteeError.getErrorInformation(PPSTEE_ERROR err);
 */

class PPSteeError {
public:

  //! Get error details
  /*!
   *
   * \param err Error code.
   * \return Human readable information on this error.
   */
  static string getErrorInformation(PPSTEE_ERROR err);

};

#endif /* PPSTEE_ERROR_HPP_ */

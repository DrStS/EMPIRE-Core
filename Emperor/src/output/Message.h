/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
/***********************************************************************************************//**
 * \file Message.h
 * This file holds the class Message
 * \date 6/8/2012
 **************************************************************************************************/

#ifndef MESSAGE_H_
#define MESSAGE_H_

#include <iostream>
#include <string>
#include <sstream>

namespace EMPIRE {
/********//**
 * \brief This manages the output functions for writing to the teminal
 **************************************************************************************************/
class Message {
public:
    /// Severity flags for Messages
    enum OutputLevel {
        ERROR,   /// Error
        WARNING, /// Warning of possible problem
        INFO,    /// Info for user
        DEBUG,   /// Debugging information in event of error
    };

    ///This variable is shared between all objects of type message
    static OutputLevel userSetOutputLevel;
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _outputLevel type enum OutputLevel set by user and stored in MetaDataBase
     * \param[in] _outputStream pointer to the ouputStream object
     * \author Stefan Sicklinger
     ***********/
    Message(OutputLevel _outputLevel, std::ostream &_outputStream);
    /***********************************************************************************************
     * \brief Destructor
     *
     * \author Stefan Sicklinger
     ***********/
    virtual ~Message();
    /***********************************************************************************************
     * \brief Operator overloading
     *
     * \author Tianyang Wang
     ***********/
    Message& operator()();
    /***********************************************************************************************
     * \brief Operator overloading
     *
     * \author Tianyang Wang
     ***********/
    Message& operator()(std::string title);
    /***********************************************************************************************
     * \brief Function template for output of basic datatypes as template parameter T
     *
     * \author Tianyang Wang
     ***********/
    template<class T> Message &operator<<(T obj) {
        if (userSetOutputLevel >= outputLevel)
            outputStream << obj;
        return (*this);
    }
    /***********************************************************************************************
     * \brief Operator overloading
     *
     * \author Tianyang Wang
     ***********/
    Message& operator<<(std::ostream& (*pf)(std::ostream&));
    /***********************************************************************************************
     * \brief Operator overloading
     *
     * \author Tianyang Wang
     ***********/
    Message& operator<<(std::ios& (*pf)(std::ios&));
    /***********************************************************************************************
     * \brief Operator overloading
     *
     * \author Tianyang Wang
     ***********/
    Message& operator<<(std::ios_base& (*pf)(std::ios_base&));
    /***********************************************************************************************
     * \brief Adapter function for Message object to std rdbuf function
     *
     * \author Stefan Sicklinger
     ***********/
    std::streambuf* rdbuf ( std::streambuf* sb ) {
        return outputStream.rdbuf(sb);
    }
    /***********************************************************************************************
     * \brief Prints heading and footing of EMPEROR output at info level
     *
     * \author Tianyang Wang, Stefan Sicklinger
     ***********/
    static void writeHeading(int headingLevel, const std::string &className,
            const std::string &info, Message &message);
    /***********************************************************************************************
     * \brief Prints text with a certain number of indents and to a certain message
     *
     * \author Tianyang Wang, Stefan Sicklinger
     ***********/
    static void writeTextWithIndent(int numIndents, const std::string &text, Message &message);
    /***********************************************************************************************
     * \brief Prints warning with heading
     *
     * \author Tianyang Wang, Stefan Sicklinger
     ***********/
    static void writeWarning(const std::string &className, const std::string &functionName,
            const std::string &message);
    /***********************************************************************************************
     * \brief Prints error with heading and calls exit()
     *
     * \author Tianyang Wang, Stefan Sicklinger
     ***********/
    static void writeError(const std::string &className, const std::string &functionName,
            const std::string &message);
    /***********************************************************************************************
     * \brief Prints EMPIRE ASCII Art
     *
     * \author Stefan Sicklinger
     ***********/
    static void writeASCIIArt();

private:
    /// outputLevel enum
    OutputLevel outputLevel;
    /// outputStream where message is redirected
    std::ostream &outputStream;
    /// outputLevelString String is appended to Message
    std::string outputLevelString;
};

extern Message infoOut;
extern Message debugOut;
extern Message errorOut;
extern Message warningOut;


/**************************************************************************************************!
      Is expanded to #pragma omp critical (IOSync);
***********/
#define CRITICAL_OUTPUT _Pragma("omp critical (IOSync)")
/**************************************************************************************************!
  Forwards \a string argument to infoOut(string) Message object;
***********/
#define INFO_OUT(string) /*
 */                CRITICAL_OUTPUT /*
 */                infoOut(string)

/**************************************************************************************************!
  Forwards \a string argument to debugOut(string) Message object;
***********/
#define DEBUG_OUT(string) /*
 */                CRITICAL_OUTPUT /*
 */                debugOut(string)
/**************************************************************************************************!
  Forwards \a string argument to errorOut(string) Message object;
***********/
#define ERROR_OUT(string) /*
 */                CRITICAL_OUTPUT /*
 */                errorOut(string)
/**************************************************************************************************!
  Forwards \a string argument to warning(string) Message object;
***********/
#define WARNING_OUT(string) /*
 */                CRITICAL_OUTPUT /*
 */                warningOut(string)
/**************************************************************************************************!
  Forwards \a headingLevel \a className \a message and \a ioObject arguments to writeHeading
  function
***********/
#define HEADING_OUT(headingLevel, className, message, ioObject) /*
*/                CRITICAL_OUTPUT /*
*/Message::writeHeading(headingLevel, className, message, ioObject)
/**************************************************************************************************!
  Forwards \a numIndents  \a message and \a ioObject arguments to writeTextWithIndent function
***********/
#define INDENT_OUT(numIndents,  message, ioObject) /*
*/                CRITICAL_OUTPUT /*
*/Message::writeTextWithIndent(numIndents, message, ioObject)
/**************************************************************************************************!
  Forwards \a className  \a functionName and \a message arguments to writeWarning function
***********/
#define WARNING_BLOCK_OUT(className,  functionName, message) /*
*/                CRITICAL_OUTPUT /*
*/Message::writeWarning(className,  functionName, message)
/**************************************************************************************************!
  Forwards \a className  \a functionName and \a message arguments to writeWarning function
***********/
#define ERROR_BLOCK_OUT(className,  functionName, message) /*
*/                CRITICAL_OUTPUT /*
*/Message::writeError(className,  functionName, message)
/**************************************************************************************************!
  Prints a ASCII Art Block(thread safe)
***********/
#define ASCIIART_BLOCK() /*
*/                CRITICAL_OUTPUT /*
*/Message::writeASCIIArt()
} /* namespace EMPIRE */
#endif /* MESSAGE_H_ */

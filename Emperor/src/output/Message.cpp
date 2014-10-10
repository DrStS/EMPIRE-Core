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
#include "Message.h"
#include <assert.h>
#include <stdlib.h>

using namespace std;

namespace EMPIRE {
///Allocate memory for static member
Message::OutputLevel Message::userSetOutputLevel = Message::DEBUG; // give a default value to it

Message debugOut  (Message::DEBUG, std::cout);
Message infoOut   (Message::INFO,  std::cout);
Message warningOut(Message::WARNING, std::cout);
Message errorOut  (Message::ERROR, std::cerr);


Message::Message(OutputLevel _outputLevel, std::ostream &_outputStream) :
        outputLevel(_outputLevel), outputStream(_outputStream) {
    if (outputLevel == ERROR)
        outputLevelString = "ERROR";
    else if (outputLevel == WARNING)
        outputLevelString = "WARNING";
    else if (outputLevel == INFO)
        outputLevelString = "INFO";
    else if (outputLevel == DEBUG)
        outputLevelString = "DEBUG";
    else
        assert(false);

    //outputStream << std::scientific;
    //outputStream.precision(10);

}

Message::~Message() {
}


Message& Message::operator()() {
    if (userSetOutputLevel >= outputLevel){
        outputStream << "EMPIRE_" << outputLevelString << ": "  ;
    }
    return (*this);
}

Message& Message::operator()(std::string title) {
    if (userSetOutputLevel >= outputLevel){
        outputStream << "EMPIRE_" << outputLevelString << ": " << title << endl;
    }
    return (*this);
}

Message& Message::operator<<(std::ostream& (*pf)(std::ostream&)) {
    if (userSetOutputLevel >= outputLevel){
        outputStream << pf;
    }
    return (*this);
}

Message& Message::operator<<(std::ios& (*pf)(std::ios&)) {
    if (userSetOutputLevel >= outputLevel){
        outputStream << pf;
    }
    return (*this);
}

Message& Message::operator<<(std::ios_base& (*pf)(std::ios_base&)) {
    if (userSetOutputLevel >= outputLevel){
        outputStream << pf;
    }
    return (*this);
}

void Message::writeHeading(int headingLevel, const std::string &className,
        const std::string &info, Message &message) {
    int numSeperationRows;
    int numCharyColumns;
    char chary;
    const int numI = 10;
    const char I = '|';
    if (headingLevel == 1) {
        numSeperationRows = 3;
        numCharyColumns = 120;
        chary = '#';
    } else if (headingLevel == 2) {
        numSeperationRows = 2;
        numCharyColumns = 120;
        chary = '=';
    } else if (headingLevel == 3) {
        numSeperationRows = 1;
        numCharyColumns = 50;
        chary = '=';
    } else if (headingLevel == 4) {
        numSeperationRows = 0;
        numCharyColumns = 0;
        chary = ' ';
    } else {
        assert(false);
    }

    for (int i = 0; i < numSeperationRows; i++)
        message << endl;
    for (int i = 0; i < numSeperationRows; i++) {
        for (int j = 0; j < numCharyColumns; j++) {
            message << chary;
        }
        message << endl;
    }
    for (int i = 0; i < numI; i++)
        message << I;
    message << " " << className << ": " << info << endl;
    for (int i = 0; i < numSeperationRows; i++) {
        for (int j = 0; j < numCharyColumns; j++) {
            message << chary;
        }
        message << endl;
    }
    for (int i = 0; i < numSeperationRows; i++)
        message << endl;
}

void Message::writeTextWithIndent(int numIndents, const string &text, Message &message) {
    string indents;
    stringstream ss(text);
    const int WORDS_PER_LINE = 10;
    for (int i = 0; i < numIndents; i++)
        indents.append("\t");

    string tmp;
    int count = 0;
    while (ss >> tmp) {
        if (count % WORDS_PER_LINE == 0) {
            if (count != 0)
                message << endl;
            message << indents;
        }
        message << tmp << " ";
        count++;
    }
    message << endl;
}


void Message::writeWarning(const string &className, const string &functionName,
        const string &message) {
    warningOut() << endl;
    warningOut() << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    warningOut() << "in " << className << "::" << functionName << "," << endl;
    warningOut() << "Message: " << message << endl;
    warningOut() << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

void Message::writeError(const string &className, const string &functionName,
        const string &message) {
    errorOut() << endl;
    errorOut() << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    errorOut() << "============================================================" << endl;
    errorOut() << "in " << className << "::" << functionName << "," << endl;
    errorOut() << "Message: " << message << endl;
    errorOut() << "============================================================" << endl;
    errorOut() << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    exit(EXIT_FAILURE);
}

void Message::writeASCIIArt() {
    infoOut  << "      .   =#+         . . .   .   . .         .   .         . .     .   . . . .       . . . . .     .     .   .     .   " << endl;
    infoOut  << "         Xx ==           .           . .   .   .           .   .   .   . .   .         .   . . .   .     . . .   .      " << endl;
    infoOut  << ". .       X+ +#+XX.     . . .   . .     .             .   . .     .   .     .             . .   .     .           . .   " << endl;
    infoOut  << "         . +x  #=xx  .             .     .     .   .   .       .       .         . .   .     .               . . .      " << endl;
    infoOut  << "   .  . +xX =X= X,#=. .   .     . .     .       . .       .           .   .           .   . . . . .   .         .     . " << endl;
    infoOut  << "  Xxxx =x. ##.  ;#+.               .   .     . . .     .   . .   . .     .       . . . .             . .     . . . .   ." << endl;
    infoOut  << "..#= +#++X= = xxx     . .     .     .     . .             . .   .         .   .     . .       . .     . .     .         " << endl;
    infoOut  << "x# xx -#  ;,+#=                .           . .         .     .     . .   .       .   .         .   .   .   .     . . .  " << endl;
    infoOut  << " +# -x- . x+-         .   .     .         . .             . . .   .     .   .           .         .   .                 " << endl;
    infoOut  << " x#X  X  XX.         .   . . .       . .         .       .         . . .   .       . . .                 . .       . .  " << endl;
    infoOut  << "  xX. . . -x=   xx#   . .   . . . .=+-+.          ;=x== . --+-=;=-=--;,   . . x=- . .=+-=;=-=-=-=,.      -+-=;=-=;+;+-=," << endl;
    infoOut  << "   -x+     .x#+#x. #=  .   . .    .#,,.# . . .    X;,x#   XX.--+-x-+-=+#+  . #+-x   =#,-=x=x-x-x-=x# .   #.-;+-+=+=+-+=#" << endl;
    infoOut  << ".    X#      =X     #x          . -x . .#   .    X+   #   #  x+=x=x+x== xx   #; # . x+ +x;+=+-+=X= ;#   ,X ;#-x++=x=x=x=" << endl;
    infoOut  << "      +#.        .+#=        xX#  .# #x X=      ,# #. x   #. #   . .  xx #;  #;.#   -# #.  .   . +X +# . # x=        .  " << endl;
    infoOut  << "    .   #x      ,Xx     . .+#+.=# .X =#- #    . # +#; #   #  X    .   -X xx  x- #   +x X- .      ,# +x  .X +x . .   .   " << endl;
    infoOut  << "         #x. .   xX,    .#Xx . .,#;# x;# x+. . x+ ++- #   #. #   .   +x. #   #; x  .=# #=-;=-=;=+x -X,   # +X=-+=+-+-x+." << endl;
    infoOut  << "        . xX= .   ;xX -x#;    ,x#,.x +,,X #   .# X +- #   #  +x=x=+=x=..#-  .x= # . +x  +-+==.,=+.Xx    ,x  ==+-+=+-=-#=" << endl;
    infoOut  << "           .x#      x#x-   .,#x+ . # #-.X.-#  x==X #,.x  .X, x-+=x=x=x+X .   #; X   =X xx+=x=#.=##   .   # +X=-=-=-=-+-." << endl;
    infoOut  << ". .     .    +#     .     +x#,    ,X.x+ ;X x  # #  X= #   # .X. , . .       .#= #   +X X-     #+ -#     ,X =+   . .     " << endl;
    infoOut  << "     . . .    ;#=  .   .x#+,   . . # X-  #- x# #;  #- x.  X, #               #, X.  =# #,. .   xX  #,    # x+    .      " << endl;
    infoOut  << ".   .         . #x    #XX   .     ,x.x+   # #+,#.  X= #   #  X               x= #   +x.x=   .   =x  #-. .x +#,=;-;=,=;=," << endl;
    infoOut  << "         .       ;x;=x         .   # #.   .X  X    X. #. .#  # . . .   .   . # .X.  .# #     .    x  #-  # . , ,., ,.- ," << endl;
}


} /* namespace EMPIRE */

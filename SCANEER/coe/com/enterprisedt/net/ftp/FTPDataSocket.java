/**
 *
 *  Java FTP client library.
 *
 *  Copyright (C) 2000-2001  Enterprise Distributed Technologies Ltd
 *
 *  www.enterprisedt.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Bug fixes, suggestions and comments should be sent to:
 *
 *  bruceb@cryptsoft.com
 *
 *  or by snail mail to:
 *
 *  Bruce P. Blackshaw
 *  53 Wakehurst Road
 *  London SW11 6DB
 *  United Kingdom
 *
 *  Change Log:
 *
 *        $Log: FTPDataSocket.java,v $
 *        Revision 1.2  2005/09/01 05:53:45  Anthony
 *        *** empty log message ***
 *
 *        Revision 1.1  2005/09/01 04:36:43  Anthony
 *        *** empty log message ***
 *
 *        Revision 1.1  2002/10/31 03:45:36  cvsuser
 *        *** empty log message ***
 *
 *        Revision 1.1  2001/10/09 20:53:46  bruceb
 *        Active mode changes
 *
 *        Revision 1.1  2001/10/05 14:42:03  bruceb
 *        moved from old project
 *
 */

package com.enterprisedt.net.ftp;

import java.io.*;
import java.net.*;

/**
 *  Supports client-side FTP DataSocket in Passive and Active Mode.
 *  Wrapper for Socket and ServerSocket. Methods are package access
 *  only - not for public use.
 *
 *  @author      Vladyslav Skarzhevskyy
 *  @version     $Revision: 1.2 $
 *
 */

public class FTPDataSocket {

    /**
     *  The underlying socket for Active connection.
     */
    private ServerSocket activeSocket = null;

    /**
     *  The underlying socket for PASV connection or Socket acepted from server.
     */
    private Socket passiveSocket = null;

    /**
     *  Create socket wrapper for Active connection.
     */
    FTPDataSocket(ServerSocket s) {
         activeSocket = s;
    }

    /**
     *  Create socket pper for PASV connection.
     */
    FTPDataSocket(Socket s) {
         passiveSocket = s;
    }


    /**
     *   Set the TCP timeout on the underlying control socket.
     *
     *   If a timeout is set, then any operation which
     *   takes longer than the timeout value will be
     *   killed with a java.io.InterruptedException.
     *
     *   @param millis The length of the timeout, in milliseconds
     */
    void setTimeout(int millis)
        throws IOException {

        if (passiveSocket != null)
            passiveSocket.setSoTimeout(millis);
        else if (activeSocket != null)
            activeSocket.setSoTimeout(millis);
    }


    /**
     *  If active mode, accepts the FTP server's connection - in PASV,
     *  we are already connected. Then gets the output stream of
     *  the connection
     *
     *  @return  output stream for underlying socket.
     */
    OutputStream getOutputStream() throws IOException {

        if (passiveSocket != null) {
            return passiveSocket.getOutputStream();
        }
        else {
            // accept socket from server, in Active mode
            passiveSocket = activeSocket.accept();
            // get and return its OutputStream
            return passiveSocket.getOutputStream ();
        }
    }

    /**
     *  If active mode, accepts the FTP server's connection - in PASV,
     *  we are already connected. Then gets the input stream of
     *  the connection
     *
     *  @return  input stream for underlying socket.
     */
    InputStream getInputStream() throws IOException {

        if (passiveSocket != null) {
            return passiveSocket.getInputStream();
        } else {
            // accept socket from server, in Active mode
            passiveSocket = activeSocket.accept();
            // get and return it's InputStream
            return passiveSocket.getInputStream ();
        }
    }

     /**
      *  Closes underlying sockets.
      */
    void close() throws IOException {

        if (passiveSocket != null)
            passiveSocket.close();
        if (activeSocket != null)
            activeSocket.close();
    }
}
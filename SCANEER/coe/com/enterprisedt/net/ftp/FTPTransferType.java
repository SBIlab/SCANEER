/**
 *  Copyright (C) 2000 Enterprise Distributed Technologies Ltd.
 *
 *
 *  Change Log:
 *
 *        $Log: FTPTransferType.java,v $
 *        Revision 1.1  2005/09/01 04:36:44  Anthony
 *        *** empty log message ***
 *
 *        Revision 1.1  2002/10/31 03:45:37  cvsuser
 *        *** empty log message ***
 *
 *        Revision 1.3  2001/10/09 20:54:08  bruceb
 *        No change
 *
 *        Revision 1.1  2001/10/05 14:42:04  bruceb
 *        moved from old project
 *
 *
 */

package com.enterprisedt.net.ftp;

/**
 *  Enumerates the transfer types possible. We
 *  support only the two common types, ASCII and
 *  Image (often called binary).
 *
 *  @author             Bruce Blackshaw
 *  @version        $Revision: 1.1 $
 *
 */
 public class FTPTransferType {

     /**
      *  Revision control id
      */
     private static String cvsId = "$Id: FTPTransferType.java,v 1.1 2005/09/01 04:36:44 Anthony Exp $";

     /**
      *   Represents ASCII transfer type
      */
     public static FTPTransferType ASCII = new FTPTransferType();

     /**
      *   Represents Image (or binary) transfer type
      */
     public static FTPTransferType BINARY = new FTPTransferType();

     /**
      *   The char sent to the server to set ASCII
      */
     static String ASCII_CHAR = "A";

     /**
      *   The char sent to the server to set BINARY
      */
     static String BINARY_CHAR = "I";

     /**
      *  Private so no-one else can instantiate this class
      */
     private FTPTransferType() {
     }
 }

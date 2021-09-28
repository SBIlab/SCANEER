/**
 *  Copyright (C) 2000 Enterprise Distributed Technologies Ltd.
 *
 *
 *  Change Log:
 *
 *        $Log: FTPConnectMode.java,v $
 *        Revision 1.1  2005/09/01 04:36:44  Anthony
 *        *** empty log message ***
 *
 *        Revision 1.1  2002/10/31 03:45:36  cvsuser
 *        *** empty log message ***
 *
 *        Revision 1.1  2001/10/09 20:53:46  bruceb
 *        Active mode changes
 *
 *        Revision 1.1  2001/10/05 14:42:04  bruceb
 *        moved from old project
 *
 *
 */

package com.enterprisedt.net.ftp;

/**
 *  Enumerates the connect modes that are possible,
 *  active & PASV
 *
 *  @author     Bruce Blackshaw
 *  @version    $Revision: 1.1 $
 *
 */
 public class FTPConnectMode {

     /**
      *  Revision control id
      */
     private static String cvsId = "$Id: FTPConnectMode.java,v 1.1 2005/09/01 04:36:44 Anthony Exp $";

     /**
      *   Represents active connect mode
      */
     public static FTPConnectMode ACTIVE = new FTPConnectMode();

     /**
      *   Represents PASV connect mode
      */
     public static FTPConnectMode PASV = new FTPConnectMode();

     /**
      *  Private so no-one else can instantiate this class
      */
     private FTPConnectMode() {
     }
 }

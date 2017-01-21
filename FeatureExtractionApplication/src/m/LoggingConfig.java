package m;

import com.sun.javafx.*;

import java.util.ArrayList;
import java.util.logging.*;
//construct a file handler and output to the mq.log file
//in the system's temp directory.



/**
 * Created by UDU on 22.07.2016.
 */

import java.io.FileInputStream;
import java.io.IOException;
import java.util.logging.*;


public class LoggingConfig {
    public static  Logger logger ;//= Logger.getLogger( LoggingConfig.class.getName());
    private static final LogManager logManager = LogManager.getLogManager();
    static{

    }

    public static void getlog(String m){

        logger  = Logger.getLogger(m);
        logger.setLevel(Level.FINEST);
        final Handler[] consoleHandler = logger.getHandlers();
        for(Handler i : consoleHandler ){
            i.setLevel(Level.FINEST);
            //i.setFormatter(new SimpleFormatter());
        }
    }



    public static void getlog(Logger logger, int level){
        if(level == 0) {
            logger.setLevel(Level.FINEST);
            ConsoleHandler consoleHandler = new ConsoleHandler();
            for(int i = 0; i < (logger.getHandlers()).length; i++){
                logger.getHandlers()[i].setLevel(Level.FINEST);
            }
        }else if(level == 5){
            logger.setLevel(Level.INFO);
//            final ConsoleHandler consoleHandler = new ConsoleHandler();
//            consoleHandler.setLevel(Level.INFO);
//            consoleHandler.setFormatter(new SimpleFormatter());
//            logger.addHandler(consoleHandler);

            ConsoleHandler consoleHandler = new ConsoleHandler();
            for(int i = 0; i < (logger.getHandlers()).length; i++){
                logger.getHandlers()[i].setLevel(Level.INFO);
            }

        }




    }


    public static void getlog_old(String m){
        logger  = Logger.getLogger(m);
        logger.setLevel(Level.FINEST);
        final ConsoleHandler consoleHandler = new ConsoleHandler();
        consoleHandler.setLevel(Level.FINEST);
        consoleHandler.setFormatter(new SimpleFormatter());
        logger.addHandler(consoleHandler);
    }

}
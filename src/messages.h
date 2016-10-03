/*
 * Messages.h
 *
 *  Created on: 29 set 2015
 *      Author: piero
 */

#ifndef MESSAGES_H_
#define MESSAGES_H_

#ifdef ENGLISH_MSG
#define MSG_WRONG_ARGC "Wrong number of input argument"
#define MSG_WRONG_ARGV "Unrecognized input argument: "
#define MSG_WARNING_ARGV "The last parameter is useless"
#define MSG_WRONG_FUNCTION "Invalid function name"
//#define MSG_WRONG_INDEX "Invalid index for decomposition counterexample function"
#else
#define MSG_WRONG_ARGC "Numero errato di argomenti di input"
#define MSG_WRONG_ARGV "Argomento di input errato: "
#define MSG_WARNING_ARGV "L'ultimo parametro e' inutile"
#define MSG_WRONG_FUNCTION "Nome della funzione non valido"
//#define MSG_WRONG_INDEX "Indice per la funzione di decomposizione del controesempio non valido" 
#endif

#endif /* MESSAGES_H_ */

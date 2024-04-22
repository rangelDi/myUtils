#include <stdio.h>
#include <stdlib.h>
#include <string.h>     // for strerror and others
#include <time.h>       // for the time stuff
#include <errno.h>      // for errno

char const chronometerFilename[] = "chronometers";
int const maxNameLength = 250;
                        
enum
{
    SUCCESS = 1,
    TOTAL_FAILURE = -1,
    TOLERABLE_FAILURE = -2
};

typedef struct
{
    struct timespec startTime;
    struct timespec endTime;
    char name[maxNameLength + 1]; // + 1 for null byte
}chronometer;

/*
 * This function returns the seconds elapsed between time two timespec structs.
 *
 * @param start_in The starting time.
 * @param end_in The ending time.
 *
 * @return The difference in seconds.
 *
 */
double secondsElapsed(struct timespec start_in, struct timespec end_in)
{
    return (end_in.tv_sec - start_in.tv_sec)
        + 1e-9 * (end_in.tv_nsec - start_in.tv_nsec);
}

/* RETURNS DYNAMICALLY ALLOCATED MEMORY
 *
 * This function attempts to read in the chronometer file, returning
 * an array of chronometer structs and the chronometer count.
 *
 * @param chronometers_out The dynamically allocated array of chronometer structs.
 * @param chronometerCount_out The number of chronometers returned.
 *
 * @return SUCCESS if a non-zero number of chronometers were read.
 *         TOTAL_FAILURE if an error occured.
 *         TOLERABLE_FAILURE if no chronometers were read because the file doesn't exist.
 */
int readChronometerFile(chronometer** chronometers_out_, int* chronometerCount_out)
{
    /* DYNAMIC DATA AND FILES */
    FILE* chronometerFile = NULL;
    chronometer* chronometers_out = NULL;

    chronometerFile = fopen(chronometerFilename,"r"); // MUST CLOSE
    if(chronometerFile == NULL)
    {
        if(errno != 2)
        {
            fprintf(stderr,"error: failed to open chronometer file for reading: %s\n", strerror(errno));
            return TOTAL_FAILURE;
        }else // error was that there was no file
        {
            return TOLERABLE_FAILURE;
        }
    }

    /* READ THE HEADER */
    // the header is a single int denoting the chronometerCount.
    int call_status = TOTAL_FAILURE;
    call_status = fread(chronometerCount_out, sizeof(int), 1, chronometerFile);
    if((*chronometerCount_out) <= 0 || call_status != 1) // there shouldn't be a 0 count
    {
        fprintf(stderr,"error: incorrect chronometer file format\n");
        fclose(chronometerFile);
        return TOTAL_FAILURE;
    }

    /* ALLOCATE MEMORY FOR READING */
    chronometers_out = malloc((*chronometerCount_out) * sizeof(chronometer)); // MUST FREE
    if(chronometers_out == NULL)
    {
        fprintf(stderr,"error: failed to allocate memory for chronometers: %s\n",strerror(errno));
        fclose(chronometerFile);
        return TOTAL_FAILURE;
    }

    /* READ CHRONOMETERS */
    call_status = fread(chronometers_out, sizeof(chronometer), (*chronometerCount_out), chronometerFile);
    if(call_status == 0 || call_status != (*chronometerCount_out))
    {
        fprintf(stderr,"error: failed to read chronometers: %s\n", strerror(errno));
        fclose(chronometerFile);
        free(chronometers_out);
        return TOTAL_FAILURE;
    }

    fclose(chronometerFile);

    (*chronometers_out_) = chronometers_out;
    return SUCCESS;
}

/*
 * This function attempts to write the array of chronometers to the file.
 * 
 * @param chronometers_in The array of chronometers to be written.
 * @param chronometerCount_in The chronometer count.
 *
 * @return SUCCESS if written successfully, TOTAL_FAILURE if any error occurs.
 */
int writeChronometerFile(chronometer* chronometers_in, int chronometerCount_in)
{
    /* DYNAMIC DATA AND FILES */
    FILE* chronometerFile = NULL;

    chronometerFile = fopen(chronometerFilename,"w"); // MUST CLOSE
    if(chronometerFile == NULL)
    {
        fprintf(stderr,"error: failed to open chronometer file for writting: %s\n", strerror(errno));
        return TOTAL_FAILURE;
    }

    /* WRITE THE HEADER */
    // the header is a single int denoting the chronometerCount.
    int call_status = TOTAL_FAILURE;
    call_status = fwrite(&chronometerCount_in, sizeof(int), 1, chronometerFile);
    if(call_status != 1) // there shouldn't be a 0 count
    {
        fprintf(stderr,"error: failed to write chronometer count: %s\n",strerror(errno));
        fclose(chronometerFile);
        return TOTAL_FAILURE;
    }

    call_status = fwrite(chronometers_in, sizeof(chronometer), chronometerCount_in, chronometerFile);
    if(call_status != chronometerCount_in)
    {
        fprintf(stderr,"error: failed to write chronometers: %s\n", strerror(errno));
        fclose(chronometerFile);
        return TOTAL_FAILURE;
    }

    fclose(chronometerFile);
    return SUCCESS;
}

/* 
 * This function will start a new chronometer and write it to file.
 *
 * @param name_in The name of the new chronometer
 * @param chronometers_in The current list of chronometers (before adding the new one)
 * @param chronometerCount_in The chronometer count before adding the new one.
 *
 * @return SUCCESS if successful or TOTAL_FAILURE on any error.
 */
int startChronometer(char* name_in, chronometer** chronometers_in_, int chronometerCount_in)
{
    chronometer* chronometers_in = (*chronometers_in_);

    if(chronometers_in != NULL) // first search to see if the chronometer exists
    {
        for(int i = 0; i < chronometerCount_in; i++)
        {
            if(strcmp(chronometers_in[i].name, name_in) == 0) // found the one we are looking for
            {
                if(chronometers_in[i].endTime.tv_sec <= 0) // if not stopped
                {
                    fprintf(stderr,"error: chronometer already running\n");
                    return TOTAL_FAILURE;
                }
                chronometers_in[i].endTime.tv_sec = -1; // set the end time to undefined again
                chronometers_in[i].endTime.tv_nsec = -1;

                // write the chronometers again
                int call_status = writeChronometerFile(chronometers_in, chronometerCount_in);
                if(call_status == TOTAL_FAILURE)
                    return TOTAL_FAILURE;

                return SUCCESS;
            }
        }
        // if we break out of the loop we didn't find it, so just continue as usual
    }

    chronometer* reallocation = realloc(chronometers_in, sizeof(chronometer)*(chronometerCount_in+1));
    if(reallocation == NULL)
    {
        fprintf(stderr,"error: failed to increment chronometer array: %s\n",strerror(errno));
        return TOTAL_FAILURE;
    }else
    {
        chronometers_in = reallocation;
        (*chronometers_in_) = chronometers_in;
    }

    int nameLength = strlen(name_in);
    if(nameLength > maxNameLength)
        nameLength = maxNameLength;

    for(int i = 0; i < nameLength; i++)
    {
        chronometers_in[chronometerCount_in].name[i] = name_in[i];
    }
    chronometers_in[chronometerCount_in].name[nameLength] = 0; // add null byte

    struct timespec startTime;
    clock_gettime(CLOCK_REALTIME, &startTime);
    chronometers_in[chronometerCount_in].startTime = startTime;

    chronometers_in[chronometerCount_in].endTime.tv_sec = -1;
    chronometers_in[chronometerCount_in].endTime.tv_nsec = -1;

    int call_status = writeChronometerFile(chronometers_in, chronometerCount_in+1); // remember we have one more chronometer now.
    if(call_status == TOTAL_FAILURE)
        return TOTAL_FAILURE;

    return SUCCESS;
}

/* 
 * This function will end a chronometer and write it to file.
 * Note that the chronometer will not be deleted from the file, it will simply be stopped.
 * Note that 'ended' chronometers may be restarted.
 *
 * @param name_in The name of the chronometer to end.
 * @param chronometers_in The current list of chronometers.
 * @param chronometerCount_in The current chronometer count.
 *
 * @return SUCCESS if successful or TOTAL_FAILURE on any error.
 */
int endChronometer(char* name_in, chronometer* chronometers_in, int chronometerCount_in)
{

    for(int i = 0; i < chronometerCount_in; i++)
    {
        if(strcmp(chronometers_in[i].name, name_in) == 0) // found the one we are looking for
        {
            struct timespec endTime;
            clock_gettime(CLOCK_REALTIME, &endTime);
            if(chronometers_in[i].endTime.tv_sec <= 0) // if there is no recorded end time
            {
                chronometers_in[i].endTime = endTime;
                int call_status = writeChronometerFile(chronometers_in, chronometerCount_in);
                if(call_status == TOTAL_FAILURE)
                    return TOTAL_FAILURE;
                else
                    return SUCCESS;
            }else
            {
                fprintf(stderr,"error: chronometer already stopped\n");
                return TOTAL_FAILURE;
            }
        }
    }
    // if we break out of this loop, we didn't find it
    fprintf(stderr,"error: no such chronometer: %s\n",name_in);
    return TOTAL_FAILURE;
}

/* 
 * This function attempts to print either all or a given chronometer in the format
 *
 *      [d]:[h]:[m]:[s]
 *
 * @param name_in The name of the chronometer to be printed or NULL to print all.
 * @param chronometers_in The array of available chronometers.
 * @param chronometerCount_in The chronometer count.
 *
 * @return TOTAL_FAILURE if a specific chronometer could not be found, SUCCESS otherwise.
 */
int printChronometer(char* name_in, chronometer* chronometers_in, int chronometerCount_in)
{
    if(name_in == NULL) //print all
    {
        for(int i = 0; i < chronometerCount_in; i++)
        {
            double elapsed;
            if(chronometers_in[i].endTime.tv_sec <= 0) // print still running
            {
                struct timespec endTime;
                clock_gettime(CLOCK_REALTIME, &endTime);
                elapsed = secondsElapsed(chronometers_in[i].startTime, endTime);
            }else // print normal
            {
                elapsed = secondsElapsed(chronometers_in[i].startTime, chronometers_in[i].endTime);
            }

            int wholeSeconds = (int)elapsed;
            double secondFraction = elapsed - wholeSeconds;

            int days = wholeSeconds/86400;
            int hours = (wholeSeconds - (days*86400))/3600;
            int minutes = (wholeSeconds - (days*86400) - (hours*3600))/60;
            float seconds = (wholeSeconds - (days*86400) - (hours*3600) - (minutes*60)) + secondFraction;

            if(chronometers_in[i].endTime.tv_sec <= 0) // print still running
            {
                printf("%s: %dd %dh %dm %fs (still running)\n",chronometers_in[i].name, days,hours,minutes,seconds);

            }else // print normal
            {
                printf("%s: %dd %dh %dm %fs (stopped)\n",chronometers_in[i].name, days,hours,minutes,seconds);
            }
        }

    }else // print just the one
    {
        for(int i = 0; i < chronometerCount_in; i++)
        {
            if(strcmp(chronometers_in[i].name, name_in) == 0) // found the one we are looking for
            {
                double elapsed;
                if(chronometers_in[i].endTime.tv_sec <= 0) // print still running
                {
                    struct timespec endTime;
                    clock_gettime(CLOCK_REALTIME, &endTime);
                    elapsed = secondsElapsed(chronometers_in[i].startTime, endTime);
                }else // print normal
                {
                    elapsed = secondsElapsed(chronometers_in[i].startTime, chronometers_in[i].endTime);
                }

                int wholeSeconds = (int)elapsed;
                double secondFraction = elapsed - wholeSeconds;

                int days = wholeSeconds/86400;
                int hours = (wholeSeconds - (days*86400))/3600;
                int minutes = (wholeSeconds - (days*86400) - (hours*3600))/60;
                float seconds = (wholeSeconds - (days*86400) - (hours*3600) - (minutes*60)) + secondFraction;

                if(chronometers_in[i].endTime.tv_sec <= 0) // print still running
                {
                    printf("%s: %dd %dh %dm %fs (still running)\n", name_in, days,hours,minutes,seconds);

                }else // print normal
                {
                    printf("%s: %dd %dh %dm %fs (stopped)\n", name_in, days,hours,minutes,seconds);
                }

                return SUCCESS;
            }
        }
        // if we break out of this loop, we didn't find it
        fprintf(stderr,"error: no such chronometer: %s\n",name_in);
        return TOTAL_FAILURE;
    }

    return SUCCESS;
}

/*
 * This function attempts to delete a chronometer from the list of chronometers and writes the file.
 * Note that the chronometerCount_in variable will not be updated because we expect program to exit after this.
 *
 * @param name_in The name of the chronometer to be deleted.
 * @param chronometers_in The array of chronometers.
 * @param chronometerCount_in The chronometerCount_in (before deletion).
 *
 * @return SUCCESS if the given chronometer is deleted, TOTAL_FAILURE otherwise.
 */
int deleteChronometer(char* name_in, chronometer* chronometers_in, int chronometerCount_in)
{
    if(name_in == NULL)
    {
        fprintf(stderr,"error: null name\n");
    }

    int target = 0;
    for(; target < chronometerCount_in; target++)
    {
        if(strcmp(chronometers_in[target].name,name_in) == 0) // found the chronometer
        {
            // shift back all the elements
            for(int i = target; i < chronometerCount_in-1; i++)
            {
                chronometers_in[i] = chronometers_in[i+1];
            }
            if(chronometerCount_in - 1 <= 0) // delete the chrono file
            {
                int call_status = remove(chronometerFilename);
                if(call_status != 0)
                {
                    fprintf(stderr, "failed to remove chronometer file: %s\n",strerror(errno));
                    return TOTAL_FAILURE;
                }else
                {
                    return SUCCESS;
                }
            }else
            {
                int call_status = writeChronometerFile(chronometers_in, chronometerCount_in-1); // write one less chronometer
                if(call_status == TOTAL_FAILURE)
                    return TOTAL_FAILURE;
                else
                    return SUCCESS;
            }

        }
    }
    // if we break out of this loop, we didn't find it
    fprintf(stderr,"error: no such chronometer: %s\n",name_in);
    return TOTAL_FAILURE;
}

/*
 * This function prints the usage of this program
 *
 * @param message An optional message to print after the usage text.
 */
void printUsage(char* message)
{
    fprintf(stderr,"usage: chrono -[option] [name]\n");
    fprintf(stderr,"option may be one of:\n");
    fprintf(stderr," Start Option: \'-s [name]\'\n");
    fprintf(stderr,"   End Option: \'-e [name]\'\n");
    fprintf(stderr," Print Option: \'-p [name]\' (this [name] is optional, if no name is provided, will print all available chronometers)\n");
    fprintf(stderr,"Delete Option: \'-d [name]\'\n");
    fprintf(stderr,"All [name] strings must be at most %d chars long\n\n", maxNameLength);

    if(message != NULL)
    {
        fprintf(stderr,"%s\n", message);
    }
}

/* 
 * This main function will take in at most one of:
 *  Start Option: '-s [name]'
 *   Stop Option: '-e [name]'
 *  Print Option: '-p [name]' (name is optional, if no name is provided, will print all available chronometers)
 * Delete Option: '-d [name]'
 */
int main(int argc, char* argv[])
{
    /* DYNAMIC DATA AND FILES */
    chronometer* chronometers = NULL;

    /* CHECK ARGS */
    if(argc > 3)
    {
        printUsage("error: too many arguments");
        exit(EXIT_FAILURE);
    }

    if(argc < 2)
    {
        printUsage("error: not enough arguments");
        exit(EXIT_FAILURE);
    }

    if((strcmp(argv[1],"-s") != 0) && (strcmp(argv[1],"-e") != 0) 
            && (strcmp(argv[1],"-p") != 0) && (strcmp(argv[1],"-d") != 0))
    {
        printUsage("error: unknown option");
        exit(EXIT_FAILURE);
    }

    /* DO SOMETHING */

    // first read the file
    int chronometerCount = 0;
    int call_status = readChronometerFile(&chronometers, &chronometerCount); // MUST FREE
    if(call_status == TOTAL_FAILURE)
    {
        exit(EXIT_FAILURE);
    }

    // print
    if(strcmp(argv[1], "-p") == 0)
    {
        if(call_status == TOLERABLE_FAILURE)
        {
            fprintf(stderr,"error: no chronometers found\n");
            exit(EXIT_FAILURE);
        }
        char* name = (argc <= 2) ? NULL : argv[2];
        call_status = printChronometer(name, chronometers, chronometerCount);

        if(call_status == SUCCESS)
        {
            free(chronometers);
            exit(EXIT_SUCCESS);
        }else
        {
            free(chronometers);
            exit(EXIT_FAILURE);
        }
    }

    if(argc <= 2)
    {
        printUsage("error: missing [name]");
        free(chronometers);
        exit(EXIT_SUCCESS);
    }

    if(strlen(argv[2]) > maxNameLength)
    {
        printUsage("error: [name] exceeds capacity");
        exit(EXIT_FAILURE);
    }

    if(strcmp(argv[1],"-s") == 0)
    {
        call_status = startChronometer(argv[2], &chronometers, chronometerCount);
    }else if(strcmp(argv[1], "-e") == 0)
    {
        call_status = endChronometer(argv[2], chronometers, chronometerCount);
    }else if(strcmp(argv[1], "-d") == 0)
    {
        call_status = deleteChronometer(argv[2], chronometers, chronometerCount);
    }else
    {
        fprintf(stderr, "????\n");
        free(chronometers);
        exit(EXIT_FAILURE);
    }

    
    if(call_status == TOTAL_FAILURE)
    {
        free(chronometers);
        exit(EXIT_FAILURE);
    }

    free(chronometers);
    return 0;
}



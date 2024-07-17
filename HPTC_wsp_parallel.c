/*************************************************************

          Parallel Branch and Bound algorithm 
      to solve the Wandering Salesman Problem (WSP)

    ------------------------------------------------------

 Author: LANGE Theo
 Date: 31/01/2023
 
 File: HPTC_wsp_parallel_s394369_Theo_LANGE.c
 Objective: This file will solve the wandering salesman problem
            using parallel a Branch and Bound algorithm
            The number of cities to be visited and the distance matrix
            will be determined using an input file.
*************************************************************/

/*******************
      Libraries
*******************/

#include <mpi.h>      // for MPI communication
#include <stdbool.h>  // for the Bool type
#include <stdio.h>    // for fprintf()
#include <stdlib.h>   // for rand() and srand()
#include <string.h>   // for memcpy()
#include <time.h>     // for time() and clock()


/*******************
  Global variables
*******************/

// Define the maximum number of cities to create array large enough for every input
#define MAX_CITIES 19

int N; // Number of cities given by the input file

int dist_matrix[MAX_CITIES][MAX_CITIES]; // Distance matrix given by the input file

int current_path[MAX_CITIES]; // Array containing the current path being studied

bool visited[MAX_CITIES]; // Boolean array, visited[i] is 0 if i has not been visited in the current path and 1 otherwise

int best_path[MAX_CITIES]; // The current best path visited

int best_path_length = 0; // The cost of the best path visited

int counter = 0; // Counter that will be used to parallelise all the different path to the processes

int new_best_path = 0; // Integer value that will be used to check wheter a better path has been found or not before synchronising the processors

double comm_t = 0.0; // Total cost of communications

/*******************
      Functions
*******************/

/*
Function: random_city
Objective: This function will randomly determine the first city of the path

Input: None
Output: return a random value between 1 and N
*/
int random_city() { return 1+rand()%N; }


/*
Function: copy_path
Objective: This function will copy the memory of an integer array of length N into a new array

Input:  Integer array dest to store the memory
        Integer array src to be copied
Output: Void
*/
void copy_path(int *dest, int *src) { memcpy(dest, src, N * sizeof(int)); }


/*
Function: path_cost
Objective: This function will copy the memory of an integer array of length N into a new array

Input:  Integer array path containing the n first visited city
        Integer n corresponding to the number of cities visited in the given path
Output: the length of the given path 
*/
int path_cost(int *path, int n) {
  int cost = 0;

  for (int i = 0; i < n - 1; i++)
    cost += dist_matrix[path[i]-1][path[i + 1]-1];

  return cost;
}


/*
Function: get_input_matrix
Objective: Create the distance matrix given an input file

Input: A pointer to the input file
Ouput: Void
*/
void get_input_matrix(FILE *file) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < i; j++) {
      fscanf(file, "%d", &dist_matrix[i][j]);
      dist_matrix[j][i] = dist_matrix[i][j];
    }
    dist_matrix[i][i] = 0;
  }
}


/*
Function: print
Objective: Print the distance matrix once created

Input: Void
Ouput: Void
*/
void print() {
  printf("Distance matrix: \n");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%d ", dist_matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}


/*
Function: min_neighbor 
Objective: Return the closest non visited city to the city given as argument

Input: A city  (Integer)
Ouput: The closest non visited city (Integer)
*/
int min_neighbor(int city) {
  int min = 0, next_city;
  // Loop over all the non visited city to determine the closest non visited one
  for (int i=0;i<N;i++) {
    if (!visited[i]) {
      if (min == 0) { // If min is null, this city is the first non visited found
        min = dist_matrix[city][i]; 
        next_city=i;
      } 
      else if (min > dist_matrix[city][i]) { 
        // update the closest city
        min = dist_matrix[city][i]; 
        next_city=i;
      } 
    }
  }

  return next_city + 1;
}


/*
Function: init_best_path 
Objective:  Initialise the best path and the bound of the problem
            The function will try to find the more efficient initial best path possible

Input: The current size of the nest path (Integer)
Ouput: Void
*/
void init_best_path(int depth) {
  if (depth != N) {
    // Find the closest non visited city to the previous city in the path
    best_path[depth] = min_neighbor(best_path[depth - 1] - 1); 
    visited[best_path[depth] - 1] = true;
    
    // Update the bound of the problem
    best_path_length += dist_matrix[best_path[depth - 1] - 1][best_path[depth] - 1];
    
    // Call itself recursively to find the next city
    init_best_path(depth+1);
    
    // reset the city to non visited
    visited[best_path[depth] - 1] = false;
  }  
}



/*
Function: branch_and_bound
Objective: Compute a recursive version Branch and Bound algorithm

Input:  The depth of the current path (Integer)
        The distance of the current path (Integer)
Ouput: Void
*/
void branch_and_bound(int depth, int current_path_length) {
  // If the current path has visited N cities
  if (depth == N) {
    if (current_path_length < best_path_length) {
      // Update the best path and its cost
      best_path_length = current_path_length;
      copy_path(best_path, current_path);
      
      // A new best path has been found
      new_best_path = 1;
    }
    return;
  } else { // else, the algorithm will look for the next city to visit
    for (int i = 0; i < N; i++) {
      if (!visited[i]) {
        visited[i] = true;
        current_path[depth] = i + 1;

        int length = current_path_length + dist_matrix[current_path[depth - 1] - 1][i];
        
        if (length < best_path_length) {
          // If the distance of the path is still lower than the bound, the algorithm will be called recursively
          branch_and_bound(depth + 1, length);
        }
        
        visited[i] = false;
      }
    }
  }
}


/*
Function: communication_point
Objective: Synchronise the best path and the length of this path between all the processes

Input:  The rank of the process (Integer)
        The number of processes used (Integer)
Ouput: Void
*/
void communication_point(int myrank, int npes) { 
  double start_t, end_t;
  
  start_t = MPI_Wtime();
  
  MPI_Allreduce(MPI_IN_PLACE, &new_best_path, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  
  // In this case, at least one process has found a better path, the processors need to be synchronised
  if (new_best_path == 1) {
    int reduce[2] = {best_path_length, myrank};
    
    // Send the minumum path cost amongst all the processes and the rank of the process thta found the new best path
    MPI_Allreduce(MPI_IN_PLACE, reduce, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
  
    // Update the bound of the problem
    best_path_length = reduce[0];
  
    // Update the best path on all the processes
    MPI_Bcast(best_path, MAX_CITIES, MPI_INT, reduce[1], MPI_COMM_WORLD);
    
    new_best_path = 0;
  }
  
  end_t = MPI_Wtime();
  comm_t += (end_t - start_t);
}


/*
Function: permute
Objective:  loop over all the available path on a certain depth to distribute them to the processes depending on their rank.
            A global counter will be used to distribute the permutations to the right process 

Input:  The depth of the current path init_depth (Integer)
        The depth that has to be reached before calling the Branch and Bound algorithm (Integer)
        The total number of initial path that has to be distributed amongst all the proceses (Integer)
        The rank of the process (Integer)
        The number of processes used (Integer)
Ouput: Void
*/
void permute(int init_depth, int depth, int nb_path, int myrank, int npes) {
    // If the depth of the current permutation is the depth to reach, the permutation is distributed to the right process
    if (init_depth == depth) {
      // If the counter has reached the number of processes , we set it to 0 again
      if (counter == npes) counter = 0;
      // When the counter is equal to the rank of the process, the process will compute the Branch and Bound method on all the route starting from the permutation
      if (counter == myrank) {
        branch_and_bound(depth, path_cost(current_path, depth));
               
        // The counter is updated
        counter += 1;
        return;
      }
      else {
        // In this case, this permutation will not be processed by this specific process, the counter is updated
        counter += 1;
        return;
      }
    }
    
    // If the depth has not been reached yet, we look for the next city to visit
    for (int j = 1; j <= N; j++) {
        if (!visited[j-1]) {
            current_path[init_depth] = j;
            visited[j-1] = true;
            permute(init_depth+1, depth, nb_path, myrank, npes);
            visited[j-1] = false;
        }
    }
}


/*
Function: main
Objective:  Compute the parallel Branch and Bound algorithm using the previous
            function

Input:  Number of command-line argument (Integer)
        String array containing the command-line argument
Ouput: Value 0
*/
int main(int argc, char *argv[]) {
  double start_t, end_t, total_t;
  int myrank, npes, result_npes, first_city, depth, nb_path, rest;
  
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Get the rank and size of the communicator
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  // Start the timer
  start_t = MPI_Wtime();

  /*********************************************
  Check whether the '-i' argument has been used
  If yes, get the name of the input file
  If no, print an error
  *********************************************/
  char *input_file = NULL;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      if (i + 1 < argc) {
        input_file = argv[i + 1];
        break;
      } else {
        printf("Error: -i option requires an argument.\n");
        return 1;
      }
    }
  }

  if (input_file) {
  } else {
    printf("Error: input file not specified.\n");
    return 1;
  }

  FILE *file;

  // Open the file and read the number of cities
  file = fopen(input_file, "r");
  fscanf(file, "%d", &N);

  // Read the distance matrix from the input file
  get_input_matrix(file);

  // Close the file
  fclose(file);

  // The process 0 will print the distance matrix and randomly pick the first_city to visit to broadcast it to all the processes
  if (myrank == 0) {
    // Print the distance matrix
    print();
    
    // Seeds the random nnumber generator used by the function rand
    srand(time(0));
    
    // Determine the first city of the route and initialise the best path and the bound of the problem
    first_city = random_city();

    
    // Broadcast the first_ciy to all the processes
    MPI_Bcast(&first_city, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    // Receive the first_city to visit from the process 0
    MPI_Bcast(&first_city, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  
  // Initialise the best path and the bound of the problem
  best_path[0] = first_city;
  visited[first_city - 1] = true;
  
  init_best_path(1);
  
  // Initialise the current path that is being processed
  current_path[0] = first_city;
  
  depth = 1;
  nb_path = N-1;

  // In order to balance the work load, the algorithm will be find the best depth to divide the tree
  while (depth != N-1 && npes > nb_path) {
    depth += 1;
    nb_path *= N-depth;
  }
  
  // Then, all the permutation at the depth determined will be distributed amongst all the processes
  if (depth == N-1) permute(1, depth, nb_path, myrank, npes);
  else {
    if (npes == 1) permute(1, depth, nb_path, myrank, npes);
    else permute(1, depth+1, nb_path, myrank, npes);
  }
  
  communication_point(myrank, npes);
  
  // End of the timer
  end_t = MPI_Wtime();
  total_t = (double)(end_t - start_t);
  
  // The process 0 prints the results
  if (myrank == 0) {
    printf("-----------------------------------------------------------------------------\n\n");
    printf("Best Route: ");
    for (int i = 0; i < N - 1; i++) {
      printf("%d -> ", best_path[i]);
    }
    printf("%d\n", best_path[N - 1]);
    printf("Tour length: %d\n\n", best_path_length);
    
    printf("-----------------------------------------------------------------------------\n\n");
    
    printf("Cost of Communications: %lf s.\n", comm_t);
    printf("Cost of Computations: %lf s.\n", total_t - comm_t);
    printf("Total Computational time: %lf s.\n\n",total_t);
  }

  MPI_Finalize();
}
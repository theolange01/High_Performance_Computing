/*************************************************************

          Serial Branch and Bound algorithm 
      to solve the Wandering Salesman Problem (WSP)

    ------------------------------------------------------

 Author: LANGE Theo
 Date: 29/01/2023
 
 File: HPTC_wsp_serial_s394369_Theo_LANGE.c
 Objective: This file will solve the wandering salesman problem
            using a serial Branch and Bound algorithm
            The number of cities to be visited and the distance matrix
            will be determined using an input file.
*************************************************************/


/*******************
      Libraries
*******************/

#include <stdbool.h> // for the Bool type
#include <stdio.h>   // for fprintf()
#include <stdlib.h>  // for rand() and srand()
#include <string.h>  // for memcpy()
#include <time.h>    // for time() and clock()


/*******************
  Global variables
*******************/

// Define the maximum number of cities to create array large enough for every input
#define MAX_CITIES 19 

int N; // Number of cities given by the input file

int dist_matrix[MAX_CITIES][MAX_CITIES]; // Distance matrix given by the input file

int current_path[MAX_CITIES]; // Array containing the current path being studied

bool visited[MAX_CITIES]; // Boolean array
// visited[i] is 0 if i has not been visited in the current path and 1 otherwise

int best_path[MAX_CITIES]; // The current best path visited

int best_path_length = 0; // The cost of the best path visited


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
Function: get_input_matrix
Objective: Create the distance matrix given an input file

Input: A pointer to the input file
Ouput: Void
*/
void get_input_matrix(FILE *file) {
  printf("Distance matrix\n");
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
      // The best path is updated if the length of the current tested route is better than the length of the best path
      best_path_length = current_path_length;
      copy_path(best_path, current_path);
    }
  } 
  else { // else, the algorithm will look for the next city to visit
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
Function: main
Objective:  Compute the serial Branch and Bound algorithm using the previous
            function

Input:  Number of command-line argument (Integer)
        String array containing the command-line argument
Ouput: Value 0
*/
int main(int argc, char *argv[]) {

  // Seeds the random nnumber generator used by the function rand
  srand(time(0));

  // Time variable
  clock_t start_t, end_t;
  double total_t;

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
      } 
      else {
        printf("Error: -i option requires an argument.\n");
        return 1;
      }
    }
  }

  if (input_file) {
  } 
  else {
    printf("Error: input file not specified.\n");
    return 1;
  }

  FILE *file;

  // Open the file and read the number of cities
  file = fopen(input_file, "r");
  fscanf(file, "%d", &N);

  // Read the distance matrix from the input file
  get_input_matrix(file);

  // Print the distance matrix
  print();

  // Close the file
  fclose(file);

  // Start the timer to determine the computational time of the serial algorithm
  start_t = clock();

  // Determine the first city of the route and initialise the best path and the bound of the problem
  int first_city = random_city();
  
  best_path[0] = first_city;
  
  current_path[0] = first_city;
  visited[first_city-1] = true;

  // Initialise the best path and the bound of the problem 
  init_best_path(1);

  // Run the branch and bound algorithm
  branch_and_bound(1, 0);

  // end the timer
  end_t = clock();

  // Determine the total computational time of the serial Branch and Bound algorithm
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;

  // Print the best tour found
  printf("Best Route: ");
  for (int i = 0; i < N - 1; i++) {
    printf("%d -> ", best_path[i]);
  }
  
  printf("%d\n\nTour length: %d\nComputational time: %lf s.\n\n", best_path[N - 1], best_path_length, total_t);
  
  return 0;
}

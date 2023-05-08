import java.io.*;
import java.util.*;

public class SAlign {
  public static void main(String[] args) {
    String input = args[0];
    String method = args[1];
    int mismatch = Integer.parseInt(args[2]);
    int gap = Integer.parseInt(args[3]);
    String out_stem = args[4];

    ArrayList<String> sequences = new ArrayList<String>(2);
    HashMap<String, ArrayList<String>> problems = new HashMap<String, ArrayList<String>>();
    String problem = "", query = "", reference = "", cigar = "";
    int[][] matrix, directions;
    int match = 0, score = 0;
    int[] scores, start_pos, end_pos;
    ArrayList<String> cigars = new ArrayList<String>();

    gap = -1 * gap;
    mismatch = -1 * mismatch;

    try {
      // Use a scanner to read the input file into a hashmap
      File inputFile = new File(input);
      Scanner sc = new Scanner(inputFile);

      while(sc.hasNextLine()) {
        problem = sc.nextLine();
        query = sc.nextLine();
        reference = sc.nextLine();
        sequences.add(query);
        sequences.add(reference);
        problems.put(problem, sequences);
        sequences = new ArrayList<String>(2);
      }

      sc.close();

    } catch (FileNotFoundException e1) {
        System.out.println("Error reading genome input file");
    }

    scores = new int[problems.size()];
    start_pos = new int[problems.size()];
    end_pos = new int[problems.size()];

    // Iterate through every problem
    for(int p = 0; p < problems.size(); p++){

      sequences = problems.get("p-" + p);

      //sequences = problems.get("p-0");
      query = sequences.get(0);
      reference = sequences.get(1);
      
      // Global sequence alignment method
      if(method.equals("global")){
        matrix = new int[query.length()+1][reference.length()+1];
        directions = new int[query.length()+1][reference.length()+1];
        directions[0][0] = 5; // 5 = done

        for(int i = 1; i <= query.length(); i++){
          matrix[i][0] = i * gap;
          directions[i][0] = 0; // 0 = gap up (in ref)
        }

        for(int i = 1; i <= reference.length(); i++){
          matrix[0][i] = i * gap;
          directions[0][i] = 1; // 1 = gap left (in query)
        }

        for(int i = 1; i <= query.length(); i++){
          for(int j = 1; j <= reference.length(); j++){
            if(query.charAt(i-1) == reference.charAt(j-1)){
              match = 0; // match cost
            } else {
              match = mismatch; // mismatch cost
            }
            
            // Minimizes cost of 3 cases
            matrix[i][j] = Math.max(Math.max(matrix[i][j-1] + gap, matrix[i-1][j] + gap), matrix[i-1][j-1] + match);
            
            // Fill in arrow matrix accordingly
            if(matrix[i][j] == matrix[i-1][j-1] + match){
              if(match == 0){
                directions[i][j] = 2; // 2 = match
              } else {
                directions[i][j] = 3; // 3 = mismatch
              }
            
            } else if(matrix[i][j] == matrix[i][j-1] + gap){
              directions[i][j] = 1; // 1 = gap left (in query)
            
            } else {
              directions[i][j] = 0; // 0 = gap up (in ref)
            }
          }
        }

        // Traceback to find cigar and score
        int i = query.length();
        int j = reference.length();
        score = 0;
        cigar = "";
        while(i > 0 || j > 0){
          // gap up
          if(directions[i][j] == 0){
            cigar = "I" + cigar;
            i = i - 1;
            score = score + gap;

          // gap left
          } else if(directions[i][j] == 1){
            cigar = "D" + cigar;
            j = j - 1;
            score = score + gap;

          // match
          } else if(directions[i][j] == 2){
            cigar = "=" + cigar;
            i = i - 1;
            j = j - 1;
          
          // mismatch
          } else if(directions[i][j] == 3){
            cigar = "X" + cigar;
            i = i - 1;
            j = j - 1;
            score = score + mismatch;
          
          // done
          } else if(directions[i][j] == 5){
            i = 0;
          }
        }

        int count = 1;
        char symbol = cigar.charAt(0);
        String new_cigar = "";

        for(int s = 1; s < cigar.length(); s++){
          if(symbol == cigar.charAt(s)){
            count++;
          } else {
            new_cigar = new_cigar + count + symbol;
            count = 1;
            symbol = cigar.charAt(s);
          }

        }
        new_cigar = new_cigar + count + symbol;
        scores[p] = score;
        cigars.add(new_cigar);

        
      
      // Local sequence alignment method
      } else {
        matrix = new int[query.length()+1][reference.length()+1];
        directions = new int[query.length()+1][reference.length()+1];
        directions[0][0] = 5; // 5 = done

        for(int i = 1; i <= query.length(); i++){
          matrix[i][0] = i * gap;
          directions[i][0] = 0; // 0 = gap up
        }

        for(int i = 1; i <= reference.length(); i++){
          matrix[0][i] = 0;
          directions[0][i] = 5; // 5 = done
        }

        for(int i = 1; i <= query.length(); i++){
          for(int j = 1; j <= reference.length(); j++){
            if(query.charAt(i-1) == reference.charAt(j-1)){
              match = 0; // match cost
            } else {
              match = mismatch; // mismatch cost
            }
            
            // Minimizes cost of 3 cases
            matrix[i][j] = Math.max(Math.max(matrix[i][j-1] + gap, matrix[i-1][j] + gap), matrix[i-1][j-1] + match);
            
            // Fill in arrow matrix accordingly
            if(matrix[i][j] == matrix[i-1][j-1] + match){
              if(match == 0){
                directions[i][j] = 2; // 2 = match
              } else {
                directions[i][j] = 3; // 3 = mismatch
              }
            
            } else if(matrix[i][j] == matrix[i][j-1] + gap){
              directions[i][j] = 1; // 1 = gap left (in query)
            
            } else {
              directions[i][j] = 0; // 0 = gap up (in ref)
            }
          }
        }

        // Find start position of query in reference
        boolean done = false;
        int j_start = 0, j_end = 0, val = -999, i_end = 0;
        //for(int i = 1; i <= query.length(); i++){
          for(int j = 1; j <= reference.length(); j++){
            
            if(matrix[query.length()][j] > val){
              val = matrix[query.length()][j];
              j_end = j;
              //i_end = i;
            }
          }
        //}

        end_pos[p] = j_end;
        System.out.println(i_end + " " + j_end);

        // Find max(opt(n,j)) to start traceback from
        
        /*for(int j = 1; j <= reference.length(); j++){
            if(matrix[query.length()][j] < val){
              val = matrix[query.length()][j];
              j_start = j;
            }
         }*/
        
        //System.out.println(j_start + " " + j_end);
        

        // Traceback to find cigar and score
        int i = query.length();
        int j = j_end;
        score = 0;
        cigar = "";
        while(i > 0 || j > 0 && matrix[i][j] != 0){
          // gap up
          if(directions[i][j] == 0){
            cigar = "I" + cigar;
            i = i - 1;
            score = score + gap;

          // gap left
          } else if(directions[i][j] == 1){
            cigar = "D" + cigar;
            j = j - 1;
            score = score + gap;

          // match
          } else if(directions[i][j] == 2){
            cigar = "=" + cigar;
            i = i - 1;
            j = j - 1;
          
          // mismatch
          } else if(directions[i][j] == 3){
            cigar = "X" + cigar;
            i = i - 1;
            j = j - 1;
            score = score + mismatch;
          
          // done
          } else if(directions[i][j] == 5){
            i = 0;
          }
        }
        start_pos[p] = j;

        int count = 1;
        char symbol = cigar.charAt(0);
        String new_cigar = "";

        for(int s = 1; s < cigar.length(); s++){
          if(symbol == cigar.charAt(s)){
            count++;
          } else {
            new_cigar = new_cigar + count + symbol;
            count = 1;
            symbol = cigar.charAt(s);
          }

        }
        new_cigar = new_cigar + count + symbol;
        scores[p] = score;
        cigars.add(new_cigar);
      }
    
    }
      
    try {

      if(method.equals("global")){
        // Construct complete output string
        StringBuilder output_string_global = new StringBuilder();

        for(int i = 0; i < problems.size(); i++){
          output_string_global.append("p-" + i + "\n");
          output_string_global.append(problems.get("p-" + i).get(0) + "\n");
          output_string_global.append(problems.get("p-" + i).get(1) + "\n");
          output_string_global.append(scores[i] + "\t0\t" + problems.get("p-" + i).get(1).length());
          output_string_global.append("\t" + cigars.get(i) + "\n");
        }

        // Create complete output file
        File out_file_global = new File(out_stem);

        if(!(out_file_global.createNewFile())) {
            System.out.println("file exists");
        }

        // Write output to complete output file
        FileWriter fw_global = new FileWriter(out_stem);
        fw_global.write(output_string_global.toString());
        fw_global.close();

      } else {
        // Construct partial output string
        StringBuilder output_string_fitting = new StringBuilder();
        for(int i = 0; i < problems.size(); i++){
          output_string_fitting.append("p-" + i + "\n");
          output_string_fitting.append(problems.get("p-" + i).get(0) + "\n");
          output_string_fitting.append(problems.get("p-" + i).get(1) + "\n");
          output_string_fitting.append(scores[i] + "\t" + start_pos[i] + "\t" + end_pos[i]);
          output_string_fitting.append("\t" + cigars.get(i) + "\n");
        }

        // Create partial output file
        File out_file_fitting = new File(out_stem);

        if(!(out_file_fitting.createNewFile())) {
            System.out.println("file exists");
        }

        // Write the output to the output file
        FileWriter fw_fitting = new FileWriter(out_stem);
        fw_fitting.write(output_string_fitting.toString());
        fw_fitting.close();
      }
      
    } catch (IOException ioe) {
      ioe.printStackTrace();
    }

  }
}

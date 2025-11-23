package com.pruebastecnicas;

import java.net.Socket;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Hello world!
 *
 */
public class App {

    /*
     * public void rotateV2(int[][] matrix) {
     * int tempSwapper = 0,tempSwapCol=0;
     * for(int i = 0;i<matrix.length;i++){
     * for(int j=0;j<matrix[i].length;j++){
     * 
     * if(j<i){
     * tempSwapper = matrix[j][i];
     * matrix[j][i] = matrix[i][j];
     * matrix[i][j] = tempSwapper;
     * 
     * }
     * }
     * 
     * 
     * 
     * }
     * 
     * for(int i = 0;i<matrix.length;i++){
     * for(int j=0;j<matrix[i].length / 2;j++){
     * int colLength = matrix[i].length-1;
     * tempSwapCol = matrix[i][j]; //[0][2],[1][2]
     * matrix[i][j] = matrix[i][colLength-j];
     * matrix[i][colLength-j] = tempSwapCol;
     * }
     * }
     * }
     */

    /**
     * public static int[] findDiagonalOrder(int[][] mat) throws
     * InterruptedException {
     * int i=0,j=0,diag=0,sizeMat = mat.length,sizeCol = mat[0].length,diagQuantity
     * = (sizeMat+sizeCol)-1,globalArr=0;
     * int[] matrixNew;
     * int kd=0;
     * if(sizeMat ==1){
     * matrixNew = new int[sizeCol];
     * 
     * while(kd<sizeCol){
     * matrixNew[kd] = mat[0][kd];
     * kd+=1;
     * }
     * return matrixNew;
     * }
     * 
     * if(sizeCol == 1){
     * matrixNew = new int[sizeMat];
     * 
     * while(kd<sizeMat){
     * matrixNew[kd] = mat[kd][0];
     * kd+=1;
     * }
     * return matrixNew;
     * }
     * 
     * matrixNew = new int[sizeMat*sizeCol];
     * //System.out.println(matrixNew.length+" tamaño arr");
     * 
     * while(diag < diagQuantity){
     * while(j<=diag){
     * int diagPairCompare = diag - j;
     * if(j<sizeMat && diagPairCompare < sizeCol){
     * matrixNew[globalArr] = mat[diagPairCompare][j];
     * globalArr++;
     * }
     * j++;
     * }
     * 
     * System.out.println("ya termino");
     * j=0;
     * diag++;
     * }
     * return matrixNew;
     * }
     * 
     */

    /*
     * public static int[][] matrixReshape(int[][] mat, int r, int c) {
     * int k=0;
     * int[][] matrixNew = new int[r][c];
     * 
     * if (mat.length * mat[0].length != r * c) return mat; // no se puede reshaping
     * for (int i = 0; i < mat.length; i++) {
     * for (int j = 0; j < mat[i].length; j++) {
     * matrixNew[k/c][k%c] = mat[i][j];
     * k++;
     * }
     * 
     * 
     * }
     * 
     * return matrixNew;
     * }
     */

    /*
     * public static boolean isToeplitzMatrix(int[][] matrix) throws
     * InterruptedException {
     * boolean finalResult = true, booleanDir = false;
     * 
     * int rowSize = matrix.length, colSize = matrix[0].length, numDiagonals =
     * (rowSize + colSize) - 1, i = 0,
     * j = colSize - 1, aux = 0, aux2 = 0, aux3 = colSize - 2;
     * 
     * while (i < numDiagonals) {
     * System.out.println(i + ">>");
     * 
     * if (!booleanDir) {
     * if (i > rowSize - 1) {
     * aux = rowSize - 1;
     * j = colSize - 2;
     * } else {
     * aux = i;
     * j = colSize - 1;
     * }
     * } else {
     * aux = rowSize - 1;
     * j = aux3;
     * }
     * 
     * while (aux > 0) {
     * if (j == colSize - 1 && i == 0 || j == 0 && i == rowSize) {
     * finalResult = true;
     * break;
     * }
     * 
     * if (j - 1 < 0) {
     * break;
     * }
     * aux2 = (aux) - 1;
     * //System.out.println(aux + ">" + j + "<" + (aux2) + "<" + (j - 1));
     * if (matrix[aux][j] == matrix[aux2][j - 1]) {
     * finalResult = true;
     * } else {
     * finalResult = false;
     * break; // aqui termina el ciclo y la validación
     * }
     * 
     * j--;
     * aux--;
     * // Thread.sleep(1000);
     * }
     * if(!finalResult) break;
     * if (i > rowSize - 1) booleanDir = true;
     * i++;
     * if (booleanDir) aux3--;
     * }
     * 
     * return finalResult;
     * }
     */

    /*
     * public static int largestOverlap(int[][] img1, int[][] img2) throws
     * InterruptedException {
     * int counter = 0, rowSize = img1.length - 1, colSize = img1[0].length - 1;
     * 
     * if (rowSize != img2.length - 1 && colSize != img2[0].length - 1)
     * return counter;
     * 
     * if (rowSize == 0 && colSize == 0 && img1[0][0] == 1 && img2[0][0] == 1) {
     * return 1;
     * } else if (rowSize == 0 && colSize == 0 && img1[0][0] == 0 && img2[0][0] ==
     * 0) {
     * return 0;
     * }
     * 
     * int x2 = 0, y2 = 0, x = 0, y = 0, i = 0, j = 0, k = 0, c = 0;
     * Set<Integer> loopStopper = new HashSet<>();
     * while (true) {
     * i = 0;
     * j = 0;
     * k = 0;
     * c = 0;
     * loopStopper.clear(); // con esto detenemos el bucle continuo
     * while (k <= img2.length - 1) {
     * while (c <= img2[0].length - 1) {
     * if (img2[k][c] == 1 && img1[k][c]!= 1) {
     * x2 = k;
     * y2 = c;
     * 
     * // System.out.println("ciclo de filas");
     * while (i <= rowSize) {
     * while (j <= colSize) {
     * x = j;
     * y = i; // al reves para match
     * if ((y2 == y ) && img1[x][y] == 1 && img2[x][y]!=1) {
     * img1[x][y] = 0;
     * img1[x2][y2] = 1;
     * counter++;
     * loopStopper.add(1);
     * } else {
     * loopStopper.add(0);
     * }
     * 
     * 
     * j += 1;
     * }
     * i += 1;
     * j = 0;
     * }
     * 
     * x = 0;
     * y = 0;
     * i = 0;
     * j = 0;
     * 
     * while (i <= rowSize) {
     * while (j <= colSize) {
     * x = i;
     * y = j; // derecho para columnas
     * if ((x2 == x ) && img1[x][y] == 1 && img1[x2][y2]!=1) {
     * img1[x][y] = 0;
     * img1[x2][y2] = 1;
     * counter++;
     * loopStopper.add(1);
     * } else {
     * loopStopper.add(0);
     * 
     * }
     * 
     * j+=1;
     * }
     * 
     * i+=1;
     * j=0;
     * }
     * 
     * i = 0;
     * j = 0;
     * while (i <= rowSize) {
     * while (j <= colSize) {
     * // si están en diagonal (diferencia absoluta de fila = de columna)
     * if (Math.abs(i - x2) == Math.abs(j - y2) && img1[i][j] == 1 && img1[x2][y2]
     * != 1) {
     * img1[i][j] = 0;
     * img1[x2][y2] = 1;
     * counter++;
     * loopStopper.add(1);
     * }
     * j++;
     * }
     * i++;
     * j = 0;
     * }
     * 
     * }else{
     * loopStopper.add(0);
     * }
     * 
     * c++;
     * }
     * k++;
     * c = 0;
     * }
     * 
     * 
     * if(!loopStopper.contains(1)){
     * break;
     * }
     * }
     * 
     * return counter;
     * }
     */

    public static int[][] transpose(int[][] matrix) {
        int i = 0, j = 0, sizeMat = matrix.length - 1, sizeCol = matrix[0].length - 1;
        int[][] visited = new int[matrix.length][matrix[0].length];

        int[][] newMatrix = new int[matrix[0].length][matrix.length];

        if (sizeMat != sizeCol) {
            for (int k = 0; k <= sizeCol; k++) {
                for (int l = 0; l <= sizeMat; l++) {

                    newMatrix[k][l] = matrix[l][k];
                }
            }
        } else {
            for (i = 0; i <= sizeMat; i++) {
                for (j = 0; j <= sizeCol; j++) {
                    if (j < sizeCol && visited[j][i] != 1) {
                        int temp = matrix[i][j];
                        matrix[i][j] = matrix[j][i];
                        matrix[j][i] = temp;

                       
                        visited[j][i] = 1;
                    }
                }
            }
        }

        return sizeMat != sizeCol ? newMatrix : matrix;
    }

    public static int longestConsecutive(int[] nums) {
        HashSet<Integer> perro = new HashSet<>((int) (nums.length / 0.75f) + 1);
        for (int num : nums)
            perro.add(num);
        // int counter=1;
        HashSet<Integer> auxBy = new HashSet<>((int) (nums.length / 0.75f) + 1);
        int numPiv = 0, longest = 0, counter = 0;
        for (int integer : nums) {
            if (!perro.contains(integer - 1))
                auxBy.add(integer);
        }

        for (Integer integer : auxBy) {
            if (!perro.contains(integer - 1)) {

                numPiv = integer;
                counter = 0;
                while (perro.contains(numPiv)) {
                    numPiv++;
                    counter++;
                }

                longest = Math.max(longest, counter);
            }
        }
        return longest;
    }

    public static String convert(String s, int numRows) throws InterruptedException {
        List<String> substrs = new ArrayList<>();

        int i = 0, j = 0;

        String subConcatter = "", finalStr = "";

        if (numRows <= 1) {
            return s;
        }

        while (i < s.length()) {
            subConcatter += s.charAt(i);
            if (subConcatter.length() == numRows) {
                substrs.add(subConcatter);
                subConcatter = "";
                if (i != s.length() - 1) {
                    subConcatter += s.charAt(i);
                }
            }

            i++;
        }

        System.out.println(i - 1);

        StringBuilder sb = new StringBuilder(numRows + subConcatter.length());
        if ((i - 1) % 2 != 0) {
            sb.append(subConcatter);
            for (int k = 0; k < numRows - 2; k++) {
                sb.append("-");
            }
        } else {
            for (int k = 0; k < numRows - 2; k++) {
                sb.append("-");
            }
            sb.append(subConcatter);
        }

        substrs.add(sb.toString());
        List<String> orderedList = new ArrayList<>();
        orderedList = IntStream
                .range(0, substrs.size())
                .mapToObj(d -> d % 2 != 0 || substrs.get(d).length() < numRows
                        ? new StringBuilder(substrs.get(d)).reverse().toString()
                        : substrs.get(d))
                .filter(h -> !(h.isBlank() && h.isEmpty()))
                .collect(Collectors.toList());
        System.out.println(orderedList.toString());

        char[] firstPivot = orderedList.get(0).toCharArray();
        for (int k = 0; k < firstPivot.length; k++) {

            System.out.println(firstPivot[k] + "PRIN");
            finalStr += firstPivot[k];

            for (int k2 = 1; k2 < orderedList.size(); k2++) {
                // if (k2%2!=0 && k<4) {
                System.out.println(orderedList.get(k2) + " index");
                if (k == firstPivot.length - 1 && (k2 % 2 == 0 || k2 % 2 != 0)) {
                    System.out.println(orderedList.get(k2).charAt(k));
                    finalStr += orderedList.get(k2).charAt(k);
                    // System.out.println(firstPivot[k]+"PRIN");
                } else {
                    if ((k == 0 || k == firstPivot.length - 1) && k2 % 2 == 0) {
                        System.out.println(orderedList.get(k2).charAt(k));
                        finalStr += orderedList.get(k2).charAt(k);
                    } else if ((k != 0 && k != firstPivot.length - 1)) {
                        System.out.println(orderedList.get(k2).charAt(k));
                        finalStr += orderedList.get(k2).charAt(k);
                    }
                }

                // }
            }
        }
        return finalStr.replace("-", "");

    }

    public static int numPairsDivisibleBy60(int[] time) {
        int counter = 0;
        int[] remainders = new int[60]; // en lugar del HashSet

        for (int t : time) {
            int rem = t % 60;
            int complement = (60 - rem) % 60; // complementario

            // aquí cuentas cuántos complementos ya existen
            counter += remainders[complement];

            // registras este residuo para futuros matches
            remainders[rem]++;
        }

        for (int i : remainders) {
            System.out.println(i);
        }

        return counter;
    }

    public static int longestPalindrome(String[] words) {
        /** */
        HashMap<String, String> cache = new HashMap<>();
        HashSet<String> pivot = new HashSet<>(Arrays.asList(words));
        HashMap<String, Integer> occurences = new HashMap();
        String build = "", middle = "";
        int count = 0, count2 = 0;
        for (int i = 0; i < words.length; i++) {
            StringBuilder sb = new StringBuilder(words[i]);
            String normal = sb.toString();
            String reverse = sb.reverse().toString();
            char[] comparison = normal.toCharArray();
            if (!cache.containsKey(normal) && pivot.contains(reverse) && (comparison[0] != comparison[1])) {
                cache.put(normal, reverse);
                cache.put(reverse, normal);
                // count+=2;
                build += normal;
            }

            if (!(comparison[0] != comparison[1])) {
                if (!occurences.containsKey(normal)) {
                    occurences.put(normal, 1);
                } else {
                    int num = occurences.get(normal);
                    num += 1;
                    occurences.put(normal, num);
                }
                // middle = normal;
            }
        }

        for (Entry<String, Integer> elem : occurences.entrySet()) {
            int resultadox = elem.getValue() / 2;

            if (resultadox > 2) {

            } else {
                middle = elem.getKey();
            }
        }

        System.out.println(occurences.toString());
        if (build.length() > 0) {
            build = build + middle + new StringBuilder(build).reverse().toString();
            /*
             * System.out.println(build);
             * System.out.println(build.equals(new
             * StringBuilder(build).reverse().toString()));
             */
            return build.length();
        } else {
            System.out.println(middle);
            return middle.length();
        }
    }

    /* primer approach 84 test cases de 102 */
    public static boolean canReorderDoubled(int[] arr) {
        boolean finalFlag = false;
        // {4,-2,2,-4}
        int i = 0;

        HashSet<Integer> searcher = new HashSet<>(Arrays.stream(arr).boxed().toList());
        List<Integer> pivot = new ArrayList<>();
        LinkedHashMap<Integer, Integer> arrFinal = new LinkedHashMap<>();
        int[] copy = Arrays.copyOf(arr, arr.length);
        Arrays.sort(arr);

        // System.out.println((arr.length-1)/2);

        while (i < arr.length / 2) {

            if ((arr[(2 * i) + 1] == 2 * arr[2 * i] || 2 * arr[(2 * i) + 1] == arr[2 * i])) {

                System.out.println("SI " + arr[(2 * i) + 1] + "<>" + 2 * arr[2 * i] + "<alm reves>"
                        + 2 * arr[(2 * i) + 1] + "<>" + arr[2 * i]);
                pivot.add(arr[2 * i]);
                pivot.add(arr[2 * i + 1]);

            }

            if (searcher.contains(arr[i] * 2)) {
                System.out.println("third case");
                arrFinal.put(arr[i] * 2, arr[i]);
            }
            i++;
        }

        i = 0;

        if (pivot.isEmpty()) {
            while (i < copy.length / 2) {
                if (copy[(2 * i) + 1] == 2 * copy[2 * i] || 2 * copy[(2 * i) + 1] == copy[2 * i]
                        || copy[(2 * i) + 1] == copy[2 * i]) {
                    System.out.println("SII " + copy[(2 * i) + 1] + "<>" + 2 * copy[2 * i] + "<alm reves>"
                            + 2 * copy[(2 * i) + 1] + "<>" + copy[2 * i]);
                    pivot.add(copy[2 * i]);
                    pivot.add(copy[2 * i + 1]);
                }

                i++;
            }
        }

        i = 0;
        if (pivot.size() == arr.length || arrFinal.size() * 2 == arr.length) {
            // System.out.println(pivot.size()+"<si>"+arr.length);
            finalFlag = true;
        }

        // System.out.println(pivot.size()+"<>"+arr.length);

        return finalFlag;
    }

    public static boolean canReorderDoubledV2(int[] arr) {
        Arrays.sort(arr);
        int i = 0;
        HashSet<Integer> tracker = new HashSet<>(Arrays.stream(arr).boxed().collect(Collectors.toList()));
        TreeMap<Integer, Integer> frequencySave = new TreeMap();
        List<Integer> newList = new ArrayList<>();
        List<Integer> errors;

        while (i < arr.length) {
            int frequency = 1;
            if (!frequencySave.containsKey(arr[i])) {
                frequencySave.put(arr[i], frequency);
            } else {
                frequency = frequencySave.get(arr[i]);
                frequency++;

                frequencySave.put(arr[i], frequency);
            }

            i++;
        }

        i = 0;

        for (Entry<Integer, Integer> mapper : frequencySave.entrySet()) {
            int c = 0;
            if (tracker.contains(mapper.getKey() * 2) && mapper.getKey() != 0) {
                int temp = Math.min(frequencySave.get(mapper.getKey()),
                        frequencySave.containsKey(mapper.getKey() * 2) ? frequencySave.get(mapper.getKey() * 2) : 0);
                while (c < temp) {
                    frequencySave.put(mapper.getKey(),
                            frequencySave.get(mapper.getKey()) != 0 ? frequencySave.get(mapper.getKey()) - 1 : 0);
                    frequencySave.put(mapper.getKey() * 2,
                            frequencySave.get(mapper.getKey() * 2) != 0 ? frequencySave.get(mapper.getKey() * 2) - 1
                                    : 0);
                    newList.add(mapper.getKey());
                    newList.add(mapper.getKey() * 2);
                    c++;
                }
            } else if (mapper.getKey() == 0 && tracker.contains(mapper.getKey() * 2) && mapper.getValue() >= 1
                    && mapper.getValue() % 2 == 0) {
                int temp = Math.min(frequencySave.get(mapper.getKey()),
                        frequencySave.containsKey(mapper.getKey() * 2) ? frequencySave.get(mapper.getKey() * 2) : 0);
                while (c < temp) {
                    frequencySave.put(mapper.getKey(),
                            frequencySave.get(mapper.getKey()) != 0 ? frequencySave.get(mapper.getKey()) - 1 : 0);
                    newList.add(mapper.getKey());
                    c++;
                }
            }
        }

        errors = frequencySave.entrySet().stream().filter(c -> c.getValue() > 0).map(y -> y.getKey())
                .collect(Collectors.toList());

        if (errors.size() > 0 || newList.size() != arr.length) {
            return false;
        } else {
            return true;
        }
    }

    public static boolean containsNearbyAlmostDuplicate(int[] nums, int indexDiff, int valueDiff) {
        // Arrays.sort(nums);

        int i = 0, j = 0;
        boolean breaker = false;
        TreeMap<Integer, Integer> finalSum = new TreeMap<>();
        while (i < nums.length) {
            while (j < nums.length) {

                if (i != j && Math.abs(i - j) <= indexDiff) {
                    // System.out.println(i+"<a>"+j+"<>"+nums[i]+"--"+nums[j]);
                    // finalSum.put(nums[j],nums[i]);
                    finalSum.put(i, j);

                }
                j++;
            }
            // System.out.println(i+"<v>"+j);

            j = 0;
            i++;
        }

        // System.out.println(finalSum.toString());

        List<Entry<Integer, Integer>> entryList = new ArrayList<>(finalSum.entrySet());
        i = finalSum.size() / 2;

        HashSet<Integer> visited = new HashSet<>();
        while (!breaker && i < finalSum.size() && i >= 0 && entryList.size() > 0) {
            // System.out.println(entryList.toString());
            Entry<Integer, Integer> mapper = entryList.get(i);
            // System.out.println(i+"]TENDRIA QUE IMPRIMIR NOMAS DOS VECES "+
            // Math.abs(nums[mapper.getKey()] - nums[mapper.getValue()] ));
            breaker = Math.abs(nums[mapper.getKey()] - nums[mapper.getValue()]) <= valueDiff ? true : false;
            // entryList.remove(i);
            // binary search

            visited.add(i);
            if ((mapper.getKey() > valueDiff || mapper.getValue() > valueDiff)) {
                i--;
            } else {
                i++;
            }
        }

        return breaker;
    }

    public static boolean isValidSudoku(char[][] board) {
        int i = 0, j = 0, k = 0, pivotter = 9, colRep = 3, rowRep = 3, pivotCol = 0,pivotRow = 0,charNumCounter=0;
        boolean returnVal = true;
        HashSet<Integer> isValid = new HashSet<>();

        HashMap<Integer,ArrayList<Map<Integer,Integer>>> coords = new HashMap<>();
        while (i < pivotter) {

            while (j < rowRep) {
                while (k < colRep) {
                    //System.out.print(board[j][k]);
                    if((int) board[j][k] != 46){
                       
                        isValid.add(Integer.valueOf(board[j][k]-'0'));
                        charNumCounter++;
                        if(!coords.containsKey(Integer.valueOf(board[j][k]-'0'))){
                         //    System.out.println("condicional primera: "+board[j][k]+" "+j+" "+k);
                            Map<Integer,Integer> coord = new HashMap();
                            coord.put(j, k);
                            
                            ArrayList<Map<Integer,Integer>> coordHand = new ArrayList<>();

                            coordHand.add(coord);
                            coords.put(Integer.valueOf( board[j][k]-'0'),coordHand);
                        }else{
                           // System.out.println("condicional alterna trigger: "+board[j][k]+" "+j+" "+k);
                            Map<Integer,Integer> coord = new HashMap();
                            coord.put(j, k);

                            ArrayList<Map<Integer,Integer>> pivot = coords.get(Integer.valueOf(board[j][k]-'0'));

                            for (Map<Integer,Integer> map : pivot) {
                                for (Entry<Integer,Integer> map2 : map.entrySet()) {
                                    if(map2.getKey() == j || map2.getValue() == k){
                                 //       System.out.println("SI DETECTA LAS COLS? " +board[j][k]);
                                        returnVal = false;
                                        break;
                                    }
                                }
                            }

                            pivot.add(coord);
                            coords.put(Integer.valueOf( board[j][k]-'0'),pivot);
                            //coords.put(Integer.valueOf(board[j][k]-'0'), .get(0).add(coord));
                        }
                    }
                    
                    k++;
                }

             //   System.out.println();
                k = pivotCol;
                j++;
            }

           // System.out.println(charNumCounter+"<>"+isValid.size());
      //     System.out.println(coords.toString());
            if(charNumCounter!=isValid.size()){
                returnVal = false;
                break;
            }

            charNumCounter=0;
            isValid.clear();

            // rowRep=0;
            colRep += 3;
            pivotCol += 3;
            k = pivotCol;
           
          //  System.out.println(i+" CUANTO PEGA ACA? " + j + " " + k + " " + rowRep + " " + colRep);
            if (i ==2 || i== 5 || i==8) {
                colRep = 3; rowRep = rowRep+=3; pivotCol = 0;pivotRow += 3;
            }
             k = pivotCol;
            j = pivotRow;
            i++;
        }

        return returnVal;
    }

    public static List<List<String>> groupAnagrams(String[] strs) {
        int i=0;
        LinkedHashMap<String,ArrayList<String>> pivot = new LinkedHashMap<>();

        List<List<String>> finalPiv = new ArrayList<>();


        for (int j = 0; j < strs.length; j++) {
            char[] piv = strs[j].toCharArray();
            char[] pivOrig = Arrays.copyOf(piv, piv.length);

            Arrays.sort(piv);

            String key = new String(piv);
            String orig = new String(pivOrig);

        //    System.out.println("iterables: "+key);
            
            if(!pivot.containsKey(key)){
                ArrayList<String> cache = new ArrayList<>();
                cache.add(orig);
                pivot.put(key, cache);
            }else{
                ArrayList<String> cache = pivot.get(key);
                cache.add(orig);
                pivot.put(key, cache);
            }
        }

        //System.out.println(pivot.toString());

        for (Entry<String,ArrayList<String>> list : pivot.entrySet()) {
            List<String> preSave = new ArrayList<>();
            for (String list2 : list.getValue()) {
         //       System.out.println(list2+" LISTA");

                preSave.add(list2);
            }

            finalPiv.add(preSave);
        }
        
     //   Collections.sort(finalPiv);
        Collections.reverse(finalPiv);
        System.out.println(finalPiv.toString());
        return finalPiv;
    }

    /*done */
    public static List<Integer> findAnagrams(String s, String p) throws InterruptedException {
        List<Integer> output = new ArrayList<Integer>();
        HashMap<Character,Integer> mapper = new HashMap<>();
        HashMap<Character,Integer> stateInit = new HashMap<>();

        for (Character letter : p.toCharArray()) {
             if(!stateInit.containsKey(letter)){
                stateInit.put(letter, 1);
            }else{
                 stateInit.put(letter, stateInit.get(letter)+1);
            }
        }
       // System.out.println(s.length()+"<>"+p.length());
        if(p.length()>s.length()){
            return output;
        }
       
        int i=0,j=0,jStop = p.length()-1;
        String nextB = s.substring(i, jStop+1);

        for (char letter : nextB.toCharArray()) {
            if(!mapper.containsKey(letter)){
                mapper.put(letter, 1);
            }else{
                 mapper.put(letter, mapper.get(letter)+1);
            }
        }

        while(jStop<s.length()){
            nextB = s.substring(i, jStop+1);
         //   System.out.println("SUBSTRING A ITERAR "+nextB);

            if(nextB.equals(p) || mapper.equals(stateInit)){
             //   System.out.println("DEBE SERVIR");
                output.add(i);
            }

            int removal = mapper.get(s.charAt(i));
            removal = removal-1 <= 0 ? 0 : removal-1;

            if(removal<=0){
                mapper.remove(s.charAt(i));
            }else{
                mapper.put(s.charAt(i), removal);
            }
            
          //  System.out.println(jStop+" valor index");
            jStop++;
            i++;
          //  System.out.println(mapper.toString()+"como queda antes");

       //     System.out.println(jStop+" valor index");
           if(jStop<s.length()){
                 mapper.put(s.charAt(jStop), mapper.containsKey(s.charAt(jStop)) ? mapper.get(s.charAt(jStop))+1 : 1);
           }

          //  System.out.println(mapper.toString()+"como queda despues");

          //  Thread.sleep(50000);
        }
        
        return output;
    }

   

/*DONE */
    public static int lengthOfLongestSubstring(String s) {
        int i = 0,last=0,j=0; boolean flag = false;

        String builder = "";

        if(s.isBlank() && !s.isEmpty()){
            return 1;
        }

        Set<String> subs = new HashSet<>();
        Set<Character> checker = new HashSet<>();
        
        while(i<s.length()){
            j=i;
            builder = "";
            while(j<s.length()){
                System.out.println(builder+"<>"+String.valueOf(s.charAt(j))+"<>"+checker.toString()+"<>"+j);
                if(!checker.contains(s.charAt(j))){
                    System.out.println("PASA");
                    builder +=String.valueOf(s.charAt(j));
                    checker.add(s.charAt(j));
                    subs.add(builder);
                }else{
                    System.out.println("no deberia "+s.charAt(j));
                    subs.add(builder);
                    flag = true;
                  //  subs.clear();
                  checker.clear();
                    builder="";
                }
                j++;

               // System.out.println(builder);
               
            }
            j=0;
            flag = false;
            i++;
        }

        subs.add(builder);

        System.out.println(subs.toString());

        System.out.println(checker.toString());

        Optional<Integer> length = subs.stream().map(charx->charx.length()).max(Comparator.naturalOrder());

        return length.isPresent() ? length.get() : 0;

    }


     public static int threeSumMulti(int[] arr, int target) throws InterruptedException {
        Arrays.sort(arr);
        TreeMap<Integer,Integer> freq = new TreeMap<>(Comparator.reverseOrder());
        for (int i : arr) {
            if (!freq.containsKey(i)) {
                freq.put(i, 1);
            }else{
                freq.put(i, freq.get(i)+1);
            }
        }

        int i=0,c=0,j=0,output=0;
        List<Integer> summer = new ArrayList<>();
        List<Integer> pivot = new ArrayList<>();

        pivot.addAll(freq.keySet());

        System.out.println(pivot.toString());
        while(i<freq.size()){
            int firstNumber = pivot.get(i);
            int firstFrequency = freq.get(firstNumber) * firstNumber;
            
            if(firstFrequency<target){
                for (int k = 0; k < freq.get(firstNumber); k++) {
                    summer.add(freq.get(firstNumber));
                }
            }else{
                summer.add(freq.get(firstNumber));
            }




            int secNum = 0;
            int thirdNum = 0;
            while(j<freq.size()){
                secNum = pivot.get(j);
                if(firstNumber+secNum < target){
                    summer.add(freq.get(pivot.get(j)));
                    /*while (c<freq.size()) {
                        thirdNum =pivot.get(c);
                        if(firstNumber+secNum+thirdNum == target){
                            summer.add(freq.get(thirdNum));
                            break;
                        }
                        c++;
                    } */

                    break;
                }     
                j++;
                c++;
            }

            i++;
          //  j=i+1;
            c=i;
            System.out.println(summer.toString());
            Thread.sleep(5000);
        }
        System.out.println(freq.toString());
        return output;
    }

    public static int[][] floodFill(int[][] image, int sr, int sc, int color)  throws InterruptedException{
        int rowSize = image.length;
        int colSize = image[0].length;
        int i=0,j=0,sr2=sr,sc2=sc;

        int firstBit = image[sr][sc];

        Set<HashMap<Integer,Integer>> visitedd = new HashSet<>();
        Set<HashMap<Integer,Integer>> coords = new HashSet<>();
        Queue<HashMap<Integer,Integer>> u1 = new LinkedList<>();
        
       
        HashMap<Integer,Integer> firstCoord = new HashMap<>(),coordCheck = new HashMap<>();
        
        firstCoord.put(sr2, sc2); // para comenzar el primer pop

        u1.add(firstCoord);
        
        while(u1.size()>0){

            if(u1.size()==0){
                break;
            }

            /**Solo en posiciones intermedias de la matriz nxm iteramos arriba, izquierda, derecha y centro, sin alterar nuestra
             * coordenada inicial para el algoritmo BFS
            */
            Iterator<Entry<Integer, Integer>> it = u1.peek().entrySet().iterator();

            while(it.hasNext()){
                Entry<Integer,Integer> valuex = it.next();

               // System.out.println(valuex.getKey()+"<>"+valuex.getValue());

                sr2 = valuex.getKey();
                sc2 = valuex.getValue();

                break;
            }
          //  sr2 = nextCoord.entry
            //arriba
            i = sr2-1;
            j = sc2;
            coordCheck.put(i, j);
            if((j<colSize && j >=0 && i < rowSize && i >=0) && !visitedd.contains(coordCheck)  && image[i][j] == firstBit){
                HashMap<Integer,Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            }else{
                HashMap<Integer,Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            //izquierda
            i = sr2;
            j = sc2-1;
            coordCheck.put(i, j);
            if((j<colSize && j >=0 && i < rowSize && i >=0) && !visitedd.contains(coordCheck) && image[i][j] == firstBit){
                HashMap<Integer,Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            }else{
                HashMap<Integer,Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            //derecha
            i = sr2;
            j = sc2+1;
            coordCheck.put(i, j);
            if((j<colSize && j >=0 && i < rowSize && i >=0) && !visitedd.contains(coordCheck) && image[i][j] == firstBit){
                HashMap<Integer,Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            }else{
                HashMap<Integer,Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            //abajo
            i = sr2+1;
            j = sc2;
            coordCheck.put(i, j);
            if((j<colSize && j >=0 && i < rowSize && i >=0) && !visitedd.contains(coordCheck) && image[i][j] == firstBit){
                HashMap<Integer,Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            }else{
                HashMap<Integer,Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            HashMap<Integer,Integer> removal = new HashMap<>();
            removal.put(sr2, sc2);
            visitedd.add(removal);

            image[sr2][sc2] = color;
      

            u1.poll(); // ahora si remueve la cabeza y ve actualizando la cola
           /* System.out.println("LA MATRIZ PERO CON EL SET DEBE TENER LAS MISMAS");
            System.out.println(visitedd.toString());
            System.out.println(coords.toString());
            System.out.println(u1.toString());
            Thread.sleep(3000); */

        }
        return image;
    }


    public static int numIslands(char[][] grid) throws InterruptedException{
        Stack<HashMap<Integer,Integer>> temporal = new Stack<>();
        Set<HashMap<Integer,Integer>> visited = new HashSet<>();

        int i=0,j=0,rowSize = grid.length,colSize=grid[0].length,x=0,y=0,counter=0;
        int[] rowDirs = {-1,0,1,0},colDirs = {0,1,0,-1};

        while(i<rowSize){

            while(j<colSize){
                if(grid[i][j]!='0'){
                 //   System.out.println("SIME DETECTA LÑOS UNOS NORMAL "+grid[i][j]);
                    HashMap<Integer,Integer> coordPivot = new HashMap<>();
                    coordPivot.put(i, j);

                    if(!visited.contains(coordPivot)){
                        temporal.add(coordPivot);
                        visited.add(coordPivot);

                        while(temporal.size()>0){
                       //     System.out.println("AQUI EMPIEZA ESTO "+temporal.toString()+"<>"+visited.toString());
                            HashMap<Integer,Integer> stackPop= temporal.pop();
                            Entry<Integer,Integer> firstCoord = stackPop.entrySet().stream().collect(Collectors.toList()).get(0);
                            visited.add(stackPop);
                          //  System.out.println();
                        //    System.out.println(firstCoord.toString());
                            for (int k = 0; k < colDirs.length; k++) {
                                x = firstCoord.getKey()+rowDirs[k];
                                y = firstCoord.getValue()+colDirs[k];

                                HashMap<Integer,Integer> coord = new HashMap<>();
                                
                          //      Thread.sleep(500);
                                if(x< rowSize && x>=0 && y<colSize && y>=0 && grid[x][y] !='0'){
                                 //   System.out.println("si lo inserta we ");
                                    coord.put(x, y);
                                    if(!visited.contains(coord)){
                                //        System.out.println("coordenada "+x+"<>"+y);
                                        temporal.add(coord);
                                    }
                                    
                                }
                            }

                        }


                     //   System.out.println("AQUI TERMINA ESTE PEDO WE "+visited.toString());

                        counter++;
                    }
                  
                }

              //  System.out.println("AQUI RETOMA EL CICLO DE NUEVO");
                j++;
            }
            j=0;
            i++;
        }

        return counter;
    }

    public static int maxAreaOfIsland(int[][] grid) throws InterruptedException{
        Stack<HashMap<Integer,Integer>> temporal = new Stack<>();
        Set<HashMap<Integer,Integer>> visited = new HashSet<>();

        int i=0,j=0,rowSize = grid.length,colSize=grid[0].length,x=0,y=0,counter=0,max=0;
        int[] rowDirs = {-1,0,1,0},colDirs = {0,1,0,-1};

        while(i<rowSize){

            while(j<colSize){
                if(grid[i][j]!=0){
                    HashMap<Integer,Integer> coordPivot = new HashMap<>();
                    coordPivot.put(i, j);

                    if(!visited.contains(coordPivot)){
                        temporal.add(coordPivot);
                        visited.add(coordPivot);
                        //cache.add('a');

                        while(temporal.size()>0){
                            HashMap<Integer,Integer> stackPop= temporal.pop();
                            
                           if(!visited.contains(stackPop)){
                                counter++;
                           }
                            Entry<Integer,Integer> firstCoord = stackPop.entrySet().stream().collect(Collectors.toList()).get(0);
                            visited.add(stackPop);
                            for (int k = 0; k < colDirs.length; k++) {
                                x = firstCoord.getKey()+rowDirs[k];
                                y = firstCoord.getValue()+colDirs[k];

                                HashMap<Integer,Integer> coord = new HashMap<>();
                                
                                if(x< rowSize && x>=0 && y<colSize && y>=0 && grid[x][y] !=0){
                                    coord.put(x, y);
                                    
                                    if(!visited.contains(coord)){
                                        temporal.add(coord);
                                    }
                                    
                                }
                            }

                        }

                        max = Math.max(max, (counter+1));
                        counter = 0;
                    }
                  
                }
                j++;
            }
            j=0;
            i++;
        }

        return max;
    }

    /*public static List<List<Integer>> subsetsWithDup(int[] nums) {
        Set<List<Integer>> output = new HashSet<>();
        int i=0,j=1;

        output.add(new ArrayList<>()); // se añade el vacio por regla
        List<Integer> mainPiv = Arrays.stream(nums).mapToObj(e->Integer.valueOf(e)).collect(Collectors.toCollection(ArrayList::new));
        Set<Integer> uniqueElements= new HashSet<>(mainPiv); //crea un set con elementos únicos
        output.add(mainPiv); //añade todo el arreglo como regla

        Iterator<Integer> it = uniqueElements.iterator();

        while(it.hasNext()){
            Integer pivot = it.next();
          //  System.out.println(pivot);
            output.add(new ArrayList<>(Arrays.asList(pivot)));
          //  output.add(new ArrayList<>(it.next()));
        }

       // System.out.println(output.toString());
        


        while(j<nums.length ){
            List<Integer> coords = new ArrayList<>();
            coords.add(nums[i]);
            coords.add(nums[j]);

            output.add(coords);
            i++;
            j++;
        }

     //   System.out.println(output.toString());

        List<List<Integer>> finalOutput = output.stream().collect(Collectors.toList());
        System.out.println(finalOutput.toString());

        

        return finalOutput;
    } */

    /* public static List<List<Integer>> subsetsWithDup(int[] nums) throws InterruptedException {
        Set<List<Integer>> output = new LinkedHashSet<>();
        int i=0,j=0;

        if(nums.length>0)output.add(new ArrayList<>()); // se añade el vacio por regla

        Queue<Integer> pivot = 

        return new ArrayList<>();
        
    } */

    public static List<List<Integer>> permute(int[] nums) throws InterruptedException{
        Queue<Integer> pivot = new LinkedList<>();
        Set<List<Integer>> output = new LinkedHashSet();

        Queue<Integer> subQueue = new LinkedList<>();
        Set<List<Integer>> subCombs = new HashSet<>();

        for (int integer : nums) {
            pivot.add(integer);
        }

       
        List<Integer> checker = pivot.stream().collect(Collectors.toList());
        List<Integer> subChequer = new ArrayList<>();
        
        while(true){
            
            output.add(new ArrayList<>(checker));
            subQueue.addAll(checker.subList(1, checker.size()));

           /*  System.out.println(subQueue.toString()+"---");
            Thread.sleep(3000); */

            subChequer.addAll(subQueue.stream().collect(Collectors.toList()));
            int localPiv = checker.get(0); 
            if(nums.length>2){
                while(!subCombs.contains(subChequer)){
                
                    Queue<Integer> cache = new LinkedList<>();
                    int rotate = subQueue.poll();
                    subQueue.add(rotate);
                    //subQueue.add(new ArrayList<>(subQueue));
                    subCombs.add(new ArrayList<>(subQueue));
                    

                    cache.add(localPiv);
                    cache.addAll(new ArrayList<>(subQueue));
                    output.add(new ArrayList<>(cache));

                    List<Integer> reversedPiv = new LinkedList<>(subQueue);
                    Collections.reverse(reversedPiv);

                    subCombs.add(reversedPiv);

                    cache.clear();
                    cache.add(localPiv);
                    cache.addAll(reversedPiv);
                    output.add(new ArrayList<>(cache));
                }
            }

            int rev = pivot.poll();
            pivot.add(rev);
            checker = pivot.stream().collect(Collectors.toList());
            if(output.contains(checker)){
                break;
            }
            subChequer.clear();
            subQueue.clear();

        }


        return output.stream().collect(Collectors.toList());
    }

    /**fibonachi hecho a la improvisada, a pesar de que ya hay un codigo, lo mismo
     * sucede con esos casos de recursión como
     */
    public static int climbStairs(int n) throws InterruptedException {
        Queue<Integer> steps = new LinkedList();
        List<List<Integer>> stepFinal = new ArrayList();
        int result = 0,counter=0,initialPiv=1;

        while(counter<n){
            System.out.println(result+"<>"+initialPiv);
          //  result=initialPiv;
           //if(result>=1){
                counter++;
            //}
            int aux = result;
            result+=initialPiv;
            if(result>=3){
                initialPiv=aux;
            }
           // initialPiv++;
       //    initialPiv=result;
            
           
        }
        System.out.println(result);
       // System.out.println(counter);
      // System.out.println(initialPiv);

        return counter;
    }

    

    public static int fib(int n) {
            if(n==0)
        {
            return 0;
        }
        else if(n==1)
        {
            return 1;
        }
        else
        {
            return fib(n-1)+fib(n-2);
        }
    }

    public static boolean canFinish(int numCourses, int[][] prerequisites) {
        Map<Integer,ArrayList<Integer>> graphStructure = new LinkedHashMap<>();
        Set<Integer> visited = new HashSet<>();
        int c=0;
        boolean result = true;
        
        if(prerequisites.length ==0) return true;

        for (int i = 0; i < prerequisites.length; i++) {
            int keyToInsert = prerequisites[i][0];
            if(!graphStructure.containsKey(keyToInsert)){
                ArrayList<Integer> pivot = new ArrayList<>();
                pivot.add(prerequisites[i][1]);
                graphStructure.put(keyToInsert, pivot);
            }else{
                ArrayList<Integer> pivot = graphStructure.get(keyToInsert);
                pivot.add(prerequisites[i][1]);
                graphStructure.put(keyToInsert, pivot);
            }
        }

        /**AVAJO */
        Iterator<Integer> keysIterator = graphStructure.keySet().iterator(); //vas a iteras los subn elementos de tu llave principal para buscarlos dentro del mismo hashmap


        System.out.println(graphStructure.toString());
        while(keysIterator.hasNext()){
            int keyToSearch = keysIterator.next(); //[[1,0],[0,1]] tienes tu grafo armado

            List<Integer> bfsSearchList = graphStructure.get(keyToSearch); // obtienes el [0] y buscas

            if(bfsSearchList.contains(keyToSearch)){ //compruebas los casos de 5=[5]
                result = false; // aunque recorras todo el grafo, ya tienes un false aunque de match en el numero de cursos con el set
            }

            while(c<bfsSearchList.size()){
                int vertexVisiting = bfsSearchList.get(c);  //1=[0,2,5,8...n] vas iterando el sub arreglo del hashmap
                if(visited.contains(vertexVisiting)){ // 2 a 1, 1 a 0 y 0 a 2 que es el keyToSearch
              
                    /*en este escenario, comprobar si al regresar está el direccionado */


                    System.out.println(keyToSearch+","+vertexVisiting+"<<<<"+visited.toString());

                    int recursivePivotKey = vertexVisiting;
                    while(visited.contains(recursivePivotKey)){ // el visited en esta fase = //[0,1]
                        Queue<Integer> greedy = new LinkedList<>(graphStructure.get(recursivePivotKey)); 

                    }
                   

                    result = false;
                    
                } //comperuebas los nodos visitados para detectar ciclos
                
                List<Integer> adjancList = graphStructure.get(vertexVisiting);

                if(adjancList!=null){
                    if(adjancList.contains(keyToSearch)){
                        result = false;
                    }
                }else{
                    visited.add(vertexVisiting);
                }
                
                c++;
            }

            visited.add(keyToSearch);

         //   System.out.println(visited.toString()+"<>"+result+"<>");

            c=0;
        }

        return visited.size() <= numCourses && result ? true:false;
    }
    

    public static void main(String[] args) throws InterruptedException {

        //{{0,10},{3,18},{5,5},{6,11},{11,14},{13,1},{15,1},{17,4}};
    //    int[][] graphs = {{1,4},{2,4},{3,1},{3,2}};

      int[][] graphs = {{1,0},{2,1},{3,2},{1,3}};
        System.out.println(canFinish(4, graphs));

      /*   System.out.println("Hello World!");
        System.out.println("ESCALERA reves "+fib(3)); */

        /*
         * int[][] prueba = { {1,2,3,4},{5,6,7,8}qr ,{9,10,11,12}};
         * 
         * System.out.println(spiralOrder(prueba));
         */
        // int[][] prueba = {{1,2,3},{4,5,6},{7,8,9}};
        // int[][] prueba = {{2,5},{8,4},{0,-1}};
        // int [] res = findDiagonalOrder(prueba);
        // int[][] prueba = { { 1, 2, 3, 4 }, { 5, 1, 2, 3 }, { 9, 5, 1, 2 } };
        // int[][] prueba = { { 1, 2 }, { 2,2 } };
        // System.out.println(isToeplitzMatrix(prueba));

        // int[][] res= matrixReshape(prueba, 1, 8);
        /**
         * for (int i = 0; i < res.length; i++) {
         * for (int j = 0; j < res[i].length; j++) {
         * System.out.println(res[i][j]);
         * }
         * }
         */
        // rotate(prueba );

        /*
         * int[][] img1 = { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
         * int[][] img2 = { { 0, 0, 0 }, { 0, 0, 0 }, { 1, 0, 0 } };
         */
        // System.out.println(largestOverlap(img1, img2));

        /**
         * int[][] img2 = {{1,2,3},{4,5,6}};
         * int [][] res= transpose(img2);
         * for (int i = 0; i < res.length; i++) {
         * for (int j = 0; j < res[i].length; j++) {
         * System.out.print(res[i][j]);
         * }
         * System.out.println();
         * }
         */

        /*
         * int[] prueba = {-1,0};
         * System.out.println(longestConsecutive(prueba));
         */

    //    String prueba = "ABC";
    //    int numRows = 2;

        // System.out.println(convert(prueba, numRows));

      //  int[] time = { 20, 40 };
        // System.out.println(numPairsDivisibleBy60(time));

        // System.out.println((60-(150%60)) %60);
        // System.out.println(150%60);

        /**
         * RETOMA ESTE PROBLERMA AL PENULTIMO O AL RESOLVER OTRO PROBLEMA
         * String[] words =
         * {"dd","aa","bb","dd","aa","dd","bb","dd","aa","cc","bb","cc","dd","cc"};
         * System.out.println(longestPalindrome(words));
         */

    //    int[] caseD = { 10, 20, 40, 80 };
        // System.out.println(Math.absExact(-3 + 2) + "+");

        // System.out.println(canReorderDoubledV2(caseD));

        /**
         * int[] arr = { 7 , 1 , 5 , 3 , 6 , 4};
         * 
         * System.out.println(returnMax(arr));
         */

     //   int[] testing = { 1, 2, 3, 1 };
        // System.out.println(Math.abs(0-3)); // es menor o igual a indexDiff que es
        // tres
   /*      char[][] arr = {{'.','.','4','.','.','.','6','3','.'},
                        {'.','.','.','.','.','.','.','.','.'},
                        {'5','.','.','.','.','.','.','9','.'},
                        {'.','.','.','5','6','.','.','.','.'},
                        {'4','.','3','.','.','.','.','.','1'},
                        {'.','.','.','7','.','.','.','.','.'},
                        {'.','.','.','5','.','.','.','.','.'},
                        {'.','.','.','.','.','.','.','.','.'},
                        {'.','.','.','.','.','.','.','.','.'}}; */
        // System.out.println(containsNearbyAlmostDuplicate(testing,3,0));

      //  System.out.println(arr[0].length + "<>" + arr.length);

     //   System.out.println(isValidSudoku(arr));

      // String[] strs = {"eat","tea","tan","ate","nat","bat"};

       //System.out.println(groupAnagrams(strs));

      /* String longest = "au";
       System.out.println(lengthOfLongestSubstring(longest)); */

    //   String s2 ="abab",p2="ab";
   //    System.out.println(findAnagrams(s, p));

     //   int[] test = {1,1,2,2,3,3,4,4,5,5};

   //     System.out.println(threeSumMulti(test, 8));

      /*   long start = System.nanoTime();
        int[][] test = {{1,1,1},{1,2,2},{1,2,2}};
        int [][] outputing = floodFill(test, 1, 1, 3);
        long end = System.nanoTime();

        for (int[] is : outputing) {
            for (int is2 : is) {
                System.out.print(is2);
            }
            System.out.println();
        }

        System.out.println((double) end-start/1000000); */

 /*        char[][] islands = {
        {'1','1','1','1','0'},
        {'1','1','0','1','0'},
        {'0','0','0','0','0'},
        {'0','0','0','0','0'}}; */



     //   System.out.println(numIslands(islands));


  /*    int[][] islandsC = {{0,0,1,0,0,0,0,1,0,0,0,0,0},
                        {0,0,0,0,0,0,0,1,1,1,0,0,0},
                        {0,1,1,0,1,0,0,0,0,0,0,0,0},
                        {0,1,0,0,1,1,0,0,1,0,1,0,0},
                        {0,1,0,0,1,1,0,0,1,1,1,0,0},
                        {0,0,0,0,0,0,0,0,0,0,1,0,0},
                        {0,0,0,0,0,0,0,1,1,1,0,0,0},
                        {0,0,0,0,0,0,0,1,1,0,0,0,0}}; */

    //System.out.println(maxAreaOfIsland(islandsC));

    int[] nums = {1,2,2};
   // System.out.println(subsetsWithDup(nums));

//   System.out.println(permute(new int[]{1,2}));

//System.out.println("ESCALERA "+climbStairs(6));


    }
}

package com.pruebastecnicas;

import java.lang.annotation.Target;
import java.math.BigInteger;
import java.net.Socket;
import java.nio.channels.Pipe.SourceChannel;
import java.time.Duration;
import java.time.Instant;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
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

    public void rotateV2(int[][] matrix) {
        int tempSwapper = 0, tempSwapCol = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {

                if (j < i) {
                    tempSwapper = matrix[j][i];
                    matrix[j][i] = matrix[i][j];
                    matrix[i][j] = tempSwapper;

                }
            }

        }

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length / 2; j++) {
                int colLength = matrix[i].length - 1;
                tempSwapCol = matrix[i][j]; // [0][2],[1][2]
                matrix[i][j] = matrix[i][colLength - j];
                matrix[i][colLength - j] = tempSwapCol;
            }
        }
    }

    public static int[] findDiagonalOrder(int[][] mat) throws InterruptedException {
        int i = 0, j = 0, diag = 0, sizeMat = mat.length, sizeCol = mat[0].length,
                diagQuantity = (sizeMat + sizeCol) - 1, globalArr = 0;
        int[] matrixNew;
        int kd = 0;
        if (sizeMat == 1) {
            matrixNew = new int[sizeCol];

            while (kd < sizeCol) {
                matrixNew[kd] = mat[0][kd];
                kd += 1;
            }
            return matrixNew;
        }

        if (sizeCol == 1) {
            matrixNew = new int[sizeMat];

            while (kd < sizeMat) {
                matrixNew[kd] = mat[kd][0];
                kd += 1;
            }
            return matrixNew;
        }
        matrixNew = new int[sizeMat * sizeCol]; // System.out.println(matrixNew.length+"
                                                // tamaño
                                                // arr");
        while (diag < diagQuantity) {
            while (j <= diag) {
                int diagPairCompare = diag - j;
                if (j < sizeMat && diagPairCompare < sizeCol) {
                    matrixNew[globalArr] = mat[diagPairCompare][j];
                    globalArr++;
                }
                j++;
            }
            System.out.println("ya termino");
            j = 0;
            diag++;
        }
        return matrixNew;
    }

    public static int[][] matrixReshape(int[][] mat, int r, int c) {
        int k = 0;
        int[][] matrixNew = new int[r][c];
        if (mat.length * mat[0].length != r * c)
            return mat; // no se puede reshaping
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[i].length; j++) {
                matrixNew[k / c][k % c] = mat[i][j];
                k++;
            }
        }
        return matrixNew;
    }

    public static boolean isToeplitzMatrix(int[][] matrix) throws InterruptedException {
        boolean finalResult = true, booleanDir = false;

        int rowSize = matrix.length, colSize = matrix[0].length, numDiagonals = (rowSize + colSize) - 1, i = 0,
                j = colSize - 1, aux = 0, aux2 = 0, aux3 = colSize - 2;

        while (i < numDiagonals) {
            System.out.println(i + ">>");

            if (!booleanDir) {
                if (i > rowSize - 1) {
                    aux = rowSize - 1;
                    j = colSize - 2;
                } else {
                    aux = i;
                    j = colSize - 1;
                }
            } else {
                aux = rowSize - 1;
                j = aux3;
            }

            while (aux > 0) {
                if (j == colSize - 1 && i == 0 || j == 0 && i == rowSize) {
                    finalResult = true;
                    break;
                }

                if (j - 1 < 0) {
                    break;
                }
                aux2 = (aux) - 1;
                // System.out.println(aux + ">" + j + "<" + (aux2) + "<" + (j - 1));
                if (matrix[aux][j] == matrix[aux2][j - 1]) {
                    finalResult = true;
                } else {
                    finalResult = false;
                    break; // aqui termina el ciclo y la validación
                }

                j--;
                aux--;
                // Thread.sleep(1000);
            }
            if (!finalResult)
                break;
            if (i > rowSize - 1)
                booleanDir = true;
            i++;
            if (booleanDir)
                aux3--;
        }

        return finalResult;
    }

    public static int largestOverlap(int[][] img1, int[][] img2) throws InterruptedException {
        int counter = 0, rowSize = img1.length - 1, colSize = img1[0].length - 1;

        if (rowSize != img2.length - 1 && colSize != img2[0].length - 1)
            return counter;

        if (rowSize == 0 && colSize == 0 && img1[0][0] == 1 && img2[0][0] == 1) {
            return 1;
        } else if (rowSize == 0 && colSize == 0 && img1[0][0] == 0 && img2[0][0] == 0) {
            return 0;
        }

        int x2 = 0, y2 = 0, x = 0, y = 0, i = 0, j = 0, k = 0, c = 0;
        Set<Integer> loopStopper = new HashSet<>();
        while (true) {
            i = 0;
            j = 0;
            k = 0;
            c = 0;
            loopStopper.clear(); // con esto detenemos el bucle continuo
            while (k <= img2.length - 1) {
                while (c <= img2[0].length - 1) {
                    if (img2[k][c] == 1 && img1[k][c] != 1) {
                        x2 = k;
                        y2 = c;

                        // System.out.println("ciclo de filas");
                        while (i <= rowSize) {
                            while (j <= colSize) {
                                x = j;
                                y = i; // al reves para match
                                if ((y2 == y) && img1[x][y] == 1 && img2[x][y] != 1) {
                                    img1[x][y] = 0;
                                    img1[x2][y2] = 1;
                                    counter++;
                                    loopStopper.add(1);
                                } else {
                                    loopStopper.add(0);
                                }

                                j += 1;
                            }
                            i += 1;
                            j = 0;
                        }

                        x = 0;
                        y = 0;
                        i = 0;
                        j = 0;

                        while (i <= rowSize) {
                            while (j <= colSize) {
                                x = i;
                                y = j; // derecho para columnas
                                if ((x2 == x) && img1[x][y] == 1 && img1[x2][y2] != 1) {
                                    img1[x][y] = 0;
                                    img1[x2][y2] = 1;
                                    counter++;
                                    loopStopper.add(1);
                                } else {
                                    loopStopper.add(0);

                                }

                                j += 1;
                            }

                            i += 1;
                            j = 0;
                        }

                        i = 0;
                        j = 0;
                        while (i <= rowSize) {
                            while (j <= colSize) {
                                // si están en diagonal (diferencia absoluta de fila = de columna)
                                if (Math.abs(i - x2) == Math.abs(j - y2) && img1[i][j] == 1 && img1[x2][y2] != 1) {
                                    img1[i][j] = 0;
                                    img1[x2][y2] = 1;
                                    counter++;
                                    loopStopper.add(1);
                                }
                                j++;
                            }
                            i++;
                            j = 0;
                        }

                    } else {
                        loopStopper.add(0);
                    }

                    c++;
                }
                k++;
                c = 0;
            }

            if (!loopStopper.contains(1)) {
                break;
            }
        }

        return counter;
    }

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
        int i = 0, j = 0, k = 0, pivotter = 9, colRep = 3, rowRep = 3, pivotCol = 0, pivotRow = 0, charNumCounter = 0;
        boolean returnVal = true;
        HashSet<Integer> isValid = new HashSet<>();

        HashMap<Integer, ArrayList<Map<Integer, Integer>>> coords = new HashMap<>();
        while (i < pivotter) {

            while (j < rowRep) {
                while (k < colRep) {
                    // System.out.print(board[j][k]);
                    if ((int) board[j][k] != 46) {

                        isValid.add(Integer.valueOf(board[j][k] - '0'));
                        charNumCounter++;
                        if (!coords.containsKey(Integer.valueOf(board[j][k] - '0'))) {
                            // System.out.println("condicional primera: "+board[j][k]+" "+j+" "+k);
                            Map<Integer, Integer> coord = new HashMap();
                            coord.put(j, k);

                            ArrayList<Map<Integer, Integer>> coordHand = new ArrayList<>();

                            coordHand.add(coord);
                            coords.put(Integer.valueOf(board[j][k] - '0'), coordHand);
                        } else {
                            // System.out.println("condicional alterna trigger: "+board[j][k]+" "+j+" "+k);
                            Map<Integer, Integer> coord = new HashMap();
                            coord.put(j, k);

                            ArrayList<Map<Integer, Integer>> pivot = coords.get(Integer.valueOf(board[j][k] - '0'));

                            for (Map<Integer, Integer> map : pivot) {
                                for (Entry<Integer, Integer> map2 : map.entrySet()) {
                                    if (map2.getKey() == j || map2.getValue() == k) {
                                        // System.out.println("SI DETECTA LAS COLS? " +board[j][k]);
                                        returnVal = false;
                                        break;
                                    }
                                }
                            }

                            pivot.add(coord);
                            coords.put(Integer.valueOf(board[j][k] - '0'), pivot);
                            // coords.put(Integer.valueOf(board[j][k]-'0'), .get(0).add(coord));
                        }
                    }

                    k++;
                }

                // System.out.println();
                k = pivotCol;
                j++;
            }

            // System.out.println(charNumCounter+"<>"+isValid.size());
            // System.out.println(coords.toString());
            if (charNumCounter != isValid.size()) {
                returnVal = false;
                break;
            }

            charNumCounter = 0;
            isValid.clear();

            // rowRep=0;
            colRep += 3;
            pivotCol += 3;
            k = pivotCol;

            // System.out.println(i+" CUANTO PEGA ACA? " + j + " " + k + " " + rowRep + " "
            // + colRep);
            if (i == 2 || i == 5 || i == 8) {
                colRep = 3;
                rowRep = rowRep += 3;
                pivotCol = 0;
                pivotRow += 3;
            }
            k = pivotCol;
            j = pivotRow;
            i++;
        }

        return returnVal;
    }

    public static List<List<String>> groupAnagrams(String[] strs) {
        int i = 0;
        LinkedHashMap<String, ArrayList<String>> pivot = new LinkedHashMap<>();

        List<List<String>> finalPiv = new ArrayList<>();

        for (int j = 0; j < strs.length; j++) {
            char[] piv = strs[j].toCharArray();
            char[] pivOrig = Arrays.copyOf(piv, piv.length);

            Arrays.sort(piv);

            String key = new String(piv);
            String orig = new String(pivOrig);

            // System.out.println("iterables: "+key);

            if (!pivot.containsKey(key)) {
                ArrayList<String> cache = new ArrayList<>();
                cache.add(orig);
                pivot.put(key, cache);
            } else {
                ArrayList<String> cache = pivot.get(key);
                cache.add(orig);
                pivot.put(key, cache);
            }
        }

        // System.out.println(pivot.toString());

        for (Entry<String, ArrayList<String>> list : pivot.entrySet()) {
            List<String> preSave = new ArrayList<>();
            for (String list2 : list.getValue()) {
                // System.out.println(list2+" LISTA");

                preSave.add(list2);
            }

            finalPiv.add(preSave);
        }

        // Collections.sort(finalPiv);
        Collections.reverse(finalPiv);
        System.out.println(finalPiv.toString());
        return finalPiv;
    }

    /* done */
    public static List<Integer> findAnagrams(String s, String p) throws InterruptedException {
        List<Integer> output = new ArrayList<Integer>();
        HashMap<Character, Integer> mapper = new HashMap<>();
        HashMap<Character, Integer> stateInit = new HashMap<>();

        for (Character letter : p.toCharArray()) {
            if (!stateInit.containsKey(letter)) {
                stateInit.put(letter, 1);
            } else {
                stateInit.put(letter, stateInit.get(letter) + 1);
            }
        }
        // System.out.println(s.length()+"<>"+p.length());
        if (p.length() > s.length()) {
            return output;
        }

        int i = 0, j = 0, jStop = p.length() - 1;
        String nextB = s.substring(i, jStop + 1);

        for (char letter : nextB.toCharArray()) {
            if (!mapper.containsKey(letter)) {
                mapper.put(letter, 1);
            } else {
                mapper.put(letter, mapper.get(letter) + 1);
            }
        }

        while (jStop < s.length()) {
            nextB = s.substring(i, jStop + 1);
            // System.out.println("SUBSTRING A ITERAR "+nextB);

            if (nextB.equals(p) || mapper.equals(stateInit)) {
                // System.out.println("DEBE SERVIR");
                output.add(i);
            }

            int removal = mapper.get(s.charAt(i));
            removal = removal - 1 <= 0 ? 0 : removal - 1;

            if (removal <= 0) {
                mapper.remove(s.charAt(i));
            } else {
                mapper.put(s.charAt(i), removal);
            }

            // System.out.println(jStop+" valor index");
            jStop++;
            i++;
            // System.out.println(mapper.toString()+"como queda antes");

            // System.out.println(jStop+" valor index");
            if (jStop < s.length()) {
                mapper.put(s.charAt(jStop), mapper.containsKey(s.charAt(jStop)) ? mapper.get(s.charAt(jStop)) + 1 : 1);
            }

            // System.out.println(mapper.toString()+"como queda despues");

            // Thread.sleep(50000);
        }

        return output;
    }

    /* DONE */
    public static int lengthOfLongestSubstring(String s) {
        int i = 0, last = 0, j = 0;
        boolean flag = false;

        String builder = "";

        if (s.isBlank() && !s.isEmpty()) {
            return 1;
        }

        Set<String> subs = new HashSet<>();
        Set<Character> checker = new HashSet<>();

        while (i < s.length()) {
            j = i;
            builder = "";
            while (j < s.length()) {
                System.out.println(builder + "<>" + String.valueOf(s.charAt(j)) + "<>" + checker.toString() + "<>" + j);
                if (!checker.contains(s.charAt(j))) {
                    System.out.println("PASA");
                    builder += String.valueOf(s.charAt(j));
                    checker.add(s.charAt(j));
                    subs.add(builder);
                } else {
                    System.out.println("no deberia " + s.charAt(j));
                    subs.add(builder);
                    flag = true;
                    // subs.clear();
                    checker.clear();
                    builder = "";
                }
                j++;

                // System.out.println(builder);

            }
            j = 0;
            flag = false;
            i++;
        }

        subs.add(builder);

        System.out.println(subs.toString());

        System.out.println(checker.toString());

        Optional<Integer> length = subs.stream().map(charx -> charx.length()).max(Comparator.naturalOrder());

        return length.isPresent() ? length.get() : 0;

    }

    public static int[][] floodFill(int[][] image, int sr, int sc, int color) throws InterruptedException {
        int rowSize = image.length;
        int colSize = image[0].length;
        int i = 0, j = 0, sr2 = sr, sc2 = sc;

        int firstBit = image[sr][sc];

        Set<HashMap<Integer, Integer>> visitedd = new HashSet<>();
        Set<HashMap<Integer, Integer>> coords = new HashSet<>();
        Queue<HashMap<Integer, Integer>> u1 = new LinkedList<>();

        HashMap<Integer, Integer> firstCoord = new HashMap<>(), coordCheck = new HashMap<>();

        firstCoord.put(sr2, sc2); // para comenzar el primer pop

        u1.add(firstCoord);

        while (u1.size() > 0) {

            if (u1.size() == 0) {
                break;
            }

            /**
             * Solo en posiciones intermedias de la matriz nxm iteramos arriba, izquierda,
             * derecha y centro, sin alterar nuestra
             * coordenada inicial para el algoritmo BFS
             */
            Iterator<Entry<Integer, Integer>> it = u1.peek().entrySet().iterator();

            while (it.hasNext()) {
                Entry<Integer, Integer> valuex = it.next();

                // System.out.println(valuex.getKey()+"<>"+valuex.getValue());

                sr2 = valuex.getKey();
                sc2 = valuex.getValue();

                break;
            }
            // sr2 = nextCoord.entry
            // arriba
            i = sr2 - 1;
            j = sc2;
            coordCheck.put(i, j);
            if ((j < colSize && j >= 0 && i < rowSize && i >= 0) && !visitedd.contains(coordCheck)
                    && image[i][j] == firstBit) {
                HashMap<Integer, Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            } else {
                HashMap<Integer, Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            // izquierda
            i = sr2;
            j = sc2 - 1;
            coordCheck.put(i, j);
            if ((j < colSize && j >= 0 && i < rowSize && i >= 0) && !visitedd.contains(coordCheck)
                    && image[i][j] == firstBit) {
                HashMap<Integer, Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            } else {
                HashMap<Integer, Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            // derecha
            i = sr2;
            j = sc2 + 1;
            coordCheck.put(i, j);
            if ((j < colSize && j >= 0 && i < rowSize && i >= 0) && !visitedd.contains(coordCheck)
                    && image[i][j] == firstBit) {
                HashMap<Integer, Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            } else {
                HashMap<Integer, Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            // abajo
            i = sr2 + 1;
            j = sc2;
            coordCheck.put(i, j);
            if ((j < colSize && j >= 0 && i < rowSize && i >= 0) && !visitedd.contains(coordCheck)
                    && image[i][j] == firstBit) {
                HashMap<Integer, Integer> pivotAdd = new HashMap<>();
                pivotAdd.put(i, j);
                coords.add(pivotAdd);
                u1.add(pivotAdd);
            } else {
                HashMap<Integer, Integer> removal = new HashMap<>();
                removal.put(i, j);
                visitedd.add(removal);
            }

            coordCheck.clear();

            HashMap<Integer, Integer> removal = new HashMap<>();
            removal.put(sr2, sc2);
            visitedd.add(removal);

            image[sr2][sc2] = color;

            u1.poll(); // ahora si remueve la cabeza y ve actualizando la cola
            /*
             * System.out.println("LA MATRIZ PERO CON EL SET DEBE TENER LAS MISMAS");
             * System.out.println(visitedd.toString());
             * System.out.println(coords.toString());
             * System.out.println(u1.toString());
             * Thread.sleep(3000);
             */

        }
        return image;
    }

    public static int numIslands(char[][] grid) throws InterruptedException {
        Stack<HashMap<Integer, Integer>> temporal = new Stack<>();
        Set<HashMap<Integer, Integer>> visited = new HashSet<>();

        int i = 0, j = 0, rowSize = grid.length, colSize = grid[0].length, x = 0, y = 0, counter = 0;
        int[] rowDirs = { -1, 0, 1, 0 }, colDirs = { 0, 1, 0, -1 };

        while (i < rowSize) {

            while (j < colSize) {
                if (grid[i][j] != '0') {
                    // System.out.println("SIME DETECTA LÑOS UNOS NORMAL "+grid[i][j]);
                    HashMap<Integer, Integer> coordPivot = new HashMap<>();
                    coordPivot.put(i, j);

                    if (!visited.contains(coordPivot)) {
                        temporal.add(coordPivot);
                        visited.add(coordPivot);

                        while (temporal.size() > 0) {
                            // System.out.println("AQUI EMPIEZA ESTO
                            // "+temporal.toString()+"<>"+visited.toString());
                            HashMap<Integer, Integer> stackPop = temporal.pop();
                            Entry<Integer, Integer> firstCoord = stackPop.entrySet().stream()
                                    .collect(Collectors.toList()).get(0);
                            visited.add(stackPop);
                            // System.out.println();
                            // System.out.println(firstCoord.toString());
                            for (int k = 0; k < colDirs.length; k++) {
                                x = firstCoord.getKey() + rowDirs[k];
                                y = firstCoord.getValue() + colDirs[k];

                                HashMap<Integer, Integer> coord = new HashMap<>();

                                // Thread.sleep(500);
                                if (x < rowSize && x >= 0 && y < colSize && y >= 0 && grid[x][y] != '0') {
                                    // System.out.println("si lo inserta we ");
                                    coord.put(x, y);
                                    if (!visited.contains(coord)) {
                                        // System.out.println("coordenada "+x+"<>"+y);
                                        temporal.add(coord);
                                    }

                                }
                            }

                        }

                        // System.out.println("AQUI TERMINA ESTE PEDO WE "+visited.toString());

                        counter++;
                    }

                }

                // System.out.println("AQUI RETOMA EL CICLO DE NUEVO");
                j++;
            }
            j = 0;
            i++;
        }

        return counter;
    }

    public static int maxAreaOfIsland(int[][] grid) throws InterruptedException {
        Stack<HashMap<Integer, Integer>> temporal = new Stack<>();
        Set<HashMap<Integer, Integer>> visited = new HashSet<>();

        int i = 0, j = 0, rowSize = grid.length, colSize = grid[0].length, x = 0, y = 0, counter = 0, max = 0;
        int[] rowDirs = { -1, 0, 1, 0 }, colDirs = { 0, 1, 0, -1 };

        while (i < rowSize) {

            while (j < colSize) {
                if (grid[i][j] != 0) {
                    HashMap<Integer, Integer> coordPivot = new HashMap<>();
                    coordPivot.put(i, j);

                    if (!visited.contains(coordPivot)) {
                        temporal.add(coordPivot);
                        visited.add(coordPivot);
                        // cache.add('a');

                        while (temporal.size() > 0) {
                            HashMap<Integer, Integer> stackPop = temporal.pop();

                            if (!visited.contains(stackPop)) {
                                counter++;
                            }
                            Entry<Integer, Integer> firstCoord = stackPop.entrySet().stream()
                                    .collect(Collectors.toList()).get(0);
                            visited.add(stackPop);
                            for (int k = 0; k < colDirs.length; k++) {
                                x = firstCoord.getKey() + rowDirs[k];
                                y = firstCoord.getValue() + colDirs[k];

                                HashMap<Integer, Integer> coord = new HashMap<>();

                                if (x < rowSize && x >= 0 && y < colSize && y >= 0 && grid[x][y] != 0) {
                                    coord.put(x, y);

                                    if (!visited.contains(coord)) {
                                        temporal.add(coord);
                                    }

                                }
                            }

                        }

                        max = Math.max(max, (counter + 1));
                        counter = 0;
                    }

                }
                j++;
            }
            j = 0;
            i++;
        }

        return max;
    }

    /**
     * fibonachi hecho a la improvisada, a pesar de que ya hay un codigo, lo mismo
     * sucede con esos casos de recursión como
     */
    public static int climbStairs(int n) throws InterruptedException {
        Queue<Integer> steps = new LinkedList();
        List<List<Integer>> stepFinal = new ArrayList();
        int result = 0, counter = 0, initialPiv = 1;

        while (counter < n) {
            System.out.println(result + "<>" + initialPiv);
            // result=initialPiv;
            // if(result>=1){
            counter++;
            // }
            int aux = result;
            result += initialPiv;
            if (result >= 3) {
                initialPiv = aux;
            }
            // initialPiv++;
            // initialPiv=result;

        }
        System.out.println(result);
        // System.out.println(counter);
        // System.out.println(initialPiv);

        return counter;
    }

    public static int fib(int n) {
        if (n == 0) {
            return 0;
        } else if (n == 1) {
            return 1;
        } else {
            return fib(n - 1) + fib(n - 2);
        }
    }

    public static boolean canFinish(int numCourses, int[][] prerequisites) throws InterruptedException {
        Map<Integer, ArrayList<Integer>> graphStructure = new LinkedHashMap<>();
        Set<Integer> visited = new HashSet<>();
        int c = 0;
        boolean result = true;

        if (prerequisites.length == 0)
            return true;

        for (int i = 0; i < prerequisites.length; i++) {
            int keyToInsert = prerequisites[i][0];
            if (!graphStructure.containsKey(keyToInsert)) {
                ArrayList<Integer> pivot = new ArrayList<>();
                pivot.add(prerequisites[i][1]);
                graphStructure.put(keyToInsert, pivot);
            } else {
                ArrayList<Integer> pivot = graphStructure.get(keyToInsert);
                pivot.add(prerequisites[i][1]);
                graphStructure.put(keyToInsert, pivot);
            }
        }

        /** AVAJO */
        Iterator<Integer> keysIterator = graphStructure.keySet().iterator(); // vas a iteras los subn elementos de tu
                                                                             // llave principal para buscarlos dentro
                                                                             // del mismo hashmap

        System.out.println(graphStructure.toString());

        // Thread.sleep(10000);
        Set<Integer> noCycleCache = new HashSet<>();
        while (keysIterator.hasNext()) {
            int keyToSearch = keysIterator.next(); // [[1,0],[0,1]] tienes tu grafo armado

            List<Integer> bfsSearchList = graphStructure.get(keyToSearch); // obtienes el [0] y buscas

            if (bfsSearchList.contains(keyToSearch)) { // compruebas los casos de 5=[5]
                // System.out.println("ACA");
                result = false; // aunque recorras todo el grafo, ya tienes un false aunque de match en el
                                // numero de cursos con el set
            }

            // comperuebas los nodos visitados para detectar ciclos
            while (c < bfsSearchList.size()) {
                int vertexVisiting = bfsSearchList.get(c); // 1=[0,2,5,8...n] vas iterando el sub arreglo del hashmap
                if (visited.contains(vertexVisiting)) { // 2 a 1, 1 a 0 y 0 a 2 que es el keyToSearch
                    /* en este escenario, comprobar si al regresar está el direccionado */
                    System.out.println(keyToSearch + "," + vertexVisiting + "<<<<EN CUAL>>>");
                    int recursivePivotKey = vertexVisiting; // 1 del 2,1
                    if (graphStructure.containsKey(recursivePivotKey)) {

                        /**
                         * checar si el cache no tiene ya identificados flujos sin ciclos
                         * para no iterar hacia atrás por cada nodo o(n)^2
                         */

                        Queue<Integer> greedy = new LinkedList<>(graphStructure.get(recursivePivotKey)); // [0,3]

                        System.out.println(greedy.toString() + "antes de filtrar" + noCycleCache.toString());
                        greedy = greedy.stream().filter(elem -> !noCycleCache.contains(elem))
                                .collect(Collectors.toCollection(LinkedList::new));
                        noCycleCache.addAll(new LinkedList<>(greedy));
                        System.out.println(greedy.toString() + "==" + noCycleCache.toString());
                        // Thread.sleep(1000);

                        if (greedy.size() == 0 && noCycleCache.contains(recursivePivotKey)) {
                            result = false;
                        }
                        if (greedy.size() > 0) {
                            while (visited.contains(recursivePivotKey)) { // el visited en esta fase = //[0,1] para
                                                                          // buscar el 1
                                if (!greedy.contains(keyToSearch)) { // si el 1 no incluye el 2 del 2,1 en el 0,3
                                    if (greedy.size() >= 1) {
                                        recursivePivotKey = greedy.poll();
                                        if (graphStructure.containsKey(recursivePivotKey)
                                                && visited.contains(recursivePivotKey)) {
                                            greedy.addAll(graphStructure.get(recursivePivotKey));
                                            noCycleCache.addAll(new LinkedList<>(greedy));
                                        }
                                    } else {
                                        break;
                                    }
                                } else {
                                    result = false;

                                    System.out.println(keyToSearch + "," + vertexVisiting
                                            + "<<<<SI DETECTA EL CIRCULAR>>>" + visited.toString());
                                    break;
                                }
                            }
                        }
                        noCycleCache.add(keyToSearch);

                    }

                }

                List<Integer> adjancList = graphStructure.get(vertexVisiting);

                if (adjancList != null) {
                    if (adjancList.contains(keyToSearch)) {
                        result = false;
                        break;
                    }
                } else {
                    visited.add(vertexVisiting);
                }

                c++;
            }

            visited.add(keyToSearch);

            if (!result) {
                break;
            }

            System.out.println(visited.toString() + "<>" + result + "<>");

            c = 0;
        }

        return visited.size() <= numCourses && result ? true : false;
    }

    public static int numberOfNeighbours(char pivot, char[][] grid, HashMap<Integer, Integer> actualCoordMap, int[] x,
            int[] y) {

        /**
         * Función auxiliar para saber el número de vecinos previo a colocar una
         * coordenada
         * en la pila de DFS, esto me hace pensar que para easte ´problema capaz
         * el BFS es más adecuado. Aunque el problema en leetcode da a pie
         * que se usan ambos.
         */
        int numberOfNeighbours = 0;

        Entry<Integer, Integer> actualCoord = actualCoordMap.entrySet().iterator().next();
        for (int i = 0; i < x.length; i++) {
            int xCord = actualCoord.getKey() + x[i];
            int yCord = actualCoord.getValue() + y[i];

            if (xCord >= 0 && xCord < grid.length && yCord >= 0 && yCord < grid[0].length) {
                if (grid[xCord][yCord] == pivot) {
                    numberOfNeighbours++;
                }
            }

        }

        return numberOfNeighbours;
    }

    public static boolean containsCycle(char[][] grid) throws InterruptedException {
        /** PRIMERO IMPLEMENTAR BFS Y DE AHÍ AGARRAR UN DFS */

        Set<HashMap<Integer, Integer>> visitedCoords = new HashSet<>(); // no usarlo en testcases gigantes!!
        /**
         * ALGO MUY INTERESANTE QUE ACABO DE CHECAR ES QUE, ES MÁS RÁPIDO
         * EL ALGORITMO PONIENDO LOS VISITED SOBRE UNA MATRIZ DEL MISMO TAMAÑO
         * QUE SOBRE UN HASHMAP QUE EN TEORIA ES MÁS RAPIDO.
         * 
         * lOGICO PUESTO QUE ES DE VELOCIDAD MÁS INSTANTANEA
         * 
         * lo cual me hizo perder tiempo de debugging extra sobre algo que ya funcionaba
         * bien
         * importante profundizar sobre las velocidades en arreglos lineales e
         * instantaneos, obvio
         * algo indexado explicitamente tiene mas velocidad que un hashcode!!!!
         * 
         */

        int[][] visited = new int[grid.length][grid[0].length];
        Stack<HashMap<Integer, Integer>> dfsLookup = new Stack<>();
        // int hasherEncoder = grid.length*grid[0].length;

        Stack<Integer> dfsLookupV2 = new Stack<>();

        int[] x = { -1, 0, 0, 1 };
        int[] y = { 0, -1, 1, 0 };

        HashMap<Integer, Integer> adder = new HashMap<>();
        int[] adderv2 = new int[2];

        int i = 0, j = 0, cycleCounter = 0, h = 0, z = 0;

        int minimumCycle = 0;

        char patternLookup = '0';

        while (h < grid.length) {

            while (z < grid[0].length) {

                patternLookup = grid[h][z]; // de acuerdo a la letra identificada veremos que es lo que resulta

                adder.clear();
                adder.put(h, z);

                adderv2[0] = h;
                adderv2[1] = z;

                int numberOfNeighbours = numberOfNeighbours(patternLookup, grid, adder, x, y);
                // System.out.println(numberOfNeighbours+" cuantos tiene");

                if (/* !visitedCoords.contains(new HashMap<>(adder)) */ visited[h][z] != 1 && numberOfNeighbours >= 2) {
                    // System.out.println(adder.toString()+" que cordenada andas iterando, para la
                    // letra "+patternLookup);

                    dfsLookupV2.add(((h * grid[0].length) + z)); // añadelo por lista y es más rapidin
                    dfsLookupV2.add(((h * grid[0].length) + z));
                    // dfsLookup.add(new HashMap<>(adder));
                    // dfsLookup.add(new HashMap<>(adder));

                    visited[h][z] = 1;
                    // visitedCoords.add(new HashMap<>(adder));

                    // si la coordenada no ha sido utilizada en algún dfs previo, continua, es
                    // irrelevante comprobar por letra, es por coordenada.

                    /*
                     * Mientras la pila esté llena, implementa el dfs a partir del punto
                     * identificado durante
                     * la matriz
                     */
                    while (dfsLookupV2.size() > 0) {

                        // Entry<Integer,Integer> coordObtainerPivot =
                        // dfsLookup.pop().entrySet().iterator().next(); // obten el primer elemento de
                        // la pila, y obten la coordenada sin usar llave

                        int coordObtainerx = dfsLookupV2.pop();

                        i = coordObtainerx / grid[0].length;
                        j = coordObtainerx % grid[0].length;

                        // i=coordObtainerPivot.getKey(); j=coordObtainerPivot.getValue();
                        // adder.clear();
                        adderv2 = new int[2];
                        // adder.put(i, j); // el visited siempre del nodo que estás recorriendo al
                        // final we

                        // visitedCoords.add(new HashMap<>(adder));
                        visited[i][j] = 1;

                        for (int k = 0; k < y.length; k++) {

                            adderv2 = new int[2];

                            /** Calcula las 4 direcciones de la coordenada */
                            int xCord = i + x[k];

                            int yCord = j + y[k];

                            /** SIGNIFICA QUE ES UNA COORDENADA VALIDA */

                            if (xCord >= 0 && xCord < grid.length && yCord >= 0 && yCord < grid[0].length) {

                                if (grid[xCord][yCord] == patternLookup) { // si el vecino es igual a la lñetra que
                                                                           // estamos buscando men...
                                    // adder.put(xCord, yCord);// procesamos la coordenada como dupla para
                                    // procesarlo más facil

                                    if (visited[xCord][yCord] != 1) { // ... y si no está en el arreglo de visitados

                                        if (dfsLookupV2.contains((xCord * grid[0].length) + yCord)) {
                                            cycleCounter++;
                                        }

                                        // dfsLookup.add(new HashMap<>(adder)); // añadelo a la pila duiramte la
                                        // iteración del for we, la coordenada vecina a la pila
                                        dfsLookupV2.add((xCord * grid[0].length) + yCord);
                                        minimumCycle++;

                                    }
                                }

                            }
                        }

                        // adder.clear();
                        adderv2 = new int[2];

                        // System.out.println(dfsLookupV2.toString()+" como es el estado de la pila
                        // después de iterar");
                        // Thread.sleep(100);
                    }

                    // System.out.println(minimumCycle+" minimum cycle");

                }
                z++;
                minimumCycle = 0;
            }
            // minimumCycle = 0;
            z = 0;
            h++;
        }

        System.out.println(cycleCounter);

        return cycleCounter > 0 ? true : false;
    }

    public static int[] maxSlidingWindowv2(int[] nums, int k) {
        long start = System.currentTimeMillis();
        int[] maxResultsPerCycle = new int[nums.length - k + 1];

        int i = 0, j = 0, arrIndex = 0;

        List<Integer> sorterHelper = new ArrayList<>();
        Deque<Integer> sortedWindow = new LinkedList<>();

        if (k > nums.length) {
            return maxResultsPerCycle; // devolver vacio para este test case
        }

        for (int l = 0; l < k; l++) {
            sorterHelper.add(nums[l]);
        }

        Collections.sort(sorterHelper, Comparator.naturalOrder());
        sortedWindow.addAll(sorterHelper);
        sorterHelper.clear();

        j = k - 1; // para el calculo de la ventana
        maxResultsPerCycle[arrIndex] = sortedWindow.getLast(); // haz la primera inserción en el arreglo
        arrIndex++;

        // System.out.println(nums[i] + "<>" + nums[j] + "<==>" +
        // sortedWindow.toString());

        /** Mientras el tope de tu ventana no toque el tope del arreglo principal */
        while (j < nums.length) {

            i++;
            j++;
            /* int[] numss = { 1, 3, -1, -3, 5, 3, 6, 7 }; */
            if (j >= nums.length)
                break;
            // System.out.println(nums[i] + "<>" + nums[j] + "<==>" +
            // sortedWindow.toString());

            if (nums[j] > sortedWindow.getLast()) {
                sortedWindow.removeLast(); // remueve el bottom de la cola y reemplazalo por el número más grande
                                           // detectado
                sortedWindow.add(nums[j]); // añade el nuevo número más grande
            } else {

                sortedWindow.add(nums[j]);
                sorterHelper.addAll(sortedWindow);
                sortedWindow.clear();
                Collections.sort(sorterHelper);

                // Collections.sort(sortedWindow);
                sortedWindow.addAll(sorterHelper);
                sorterHelper.clear();
                // System.out.println("aca we se deberia quitar eso "+sortedWindow.toString());
            }

            if (nums[i] == sortedWindow.getFirst()) {
                sortedWindow.removeFirst();
            }
            maxResultsPerCycle[arrIndex] = sortedWindow.getLast(); // haz la primera inserción en el arreglo
            arrIndex++;
        }

        long end = System.currentTimeMillis();

        for (int integer : maxResultsPerCycle) {
            System.out.println(integer);
        }

        return maxResultsPerCycle;
    }

    public static int[] maxSlidingWindow(int[] nums, int k) throws InterruptedException {
        // Instant start = Instant.now();
        int[] maxResultsPerCycle = new int[nums.length - k + 1];

        int i = 0, j = 0, arrIndex = 0;

        Deque<Integer> sortedWindow = new LinkedList<>();

        if (k > nums.length) {
            return maxResultsPerCycle; // devolver vacio para este test case
        }

        /** 7,2,4 = 2 && 9,11 2 */

        if (k == 1) {
            return nums;
        }

        /**
         * Estructura el comportamiento de la cola monotona desde el for de acuerdo al
         * primer rango
         */
        for (int l = 0; l < k; l++) {

            /** Esta es la lógica de un monotonic queue */
            if (sortedWindow.size() == 0) {
                sortedWindow.add(l);
            } else {
                while (sortedWindow.size() > 0 && nums[l] >= nums[sortedWindow.getLast()]) {
                    sortedWindow.removeLast(); // comportamiento de una cola monotona
                }
                sortedWindow.addLast(l);

            }
        }

        System.out.println(sortedWindow.toString());
        // Thread.sleep(10000);
        j = k - 1;

        while (j < nums.length) {
            // System.out.println(sortedWindow.toString()+" estado de la pila
            // "+nums[i]+"-"+nums[j] +"->>"+i+""+j);
            if (i > 0) {

                while (sortedWindow.size() > 0 && nums[j] >= nums[sortedWindow.getLast()]) {
                    sortedWindow.removeLast(); // comportamiento de una cola monotona
                }

                sortedWindow.addLast(j);

            }

            maxResultsPerCycle[arrIndex] = nums[sortedWindow.getFirst()]; // haz la primera inserción en el arreglo
            arrIndex++;

            if (i == sortedWindow.getFirst()) {
                sortedWindow.removeFirst();
            }

            /** Haz crecer la ventana en todo momento */
            i++;
            j++;
        }

        /*
         * for (int integer : maxResultsPerCycle) {
         * System.out.println(integer);
         * }
         */
        return maxResultsPerCycle;
    }

    // target = 7, nums = [2,3,1,2,4,3]
    public static int minSubArrayLen(int target, int[] nums) throws InterruptedException {
        int izq = 0, der = 0, indexCounter = 0, sumAdder = 0, finalIndex = Integer.MAX_VALUE;
        HashMap<Integer, Boolean> hasCheck = new HashMap<>();

        while (der < nums.length) {
            sumAdder += nums[der];

            if (sumAdder >= target) {
                // System.out.println("INCREMENTO "+sumAdder+"<>"+izq+" "+der);
                finalIndex = Math.min(Math.abs(izq - der) + 1, finalIndex);
                hasCheck.put(target, true);

                while (sumAdder > target) {
                    sumAdder -= nums[izq];
                    izq++;
                    indexCounter--;
                    if (sumAdder >= target) {
                        finalIndex = Math.min(Math.abs(izq - der) + 1, finalIndex);
                        hasCheck.put(target, true);
                        // System.out.println("deCREMENTO "+sumAdder+"<>"+izq+" "+der);
                    }
                }
            }

            der++;
            indexCounter++;

        }

        return hasCheck.containsKey(target) ? finalIndex : 0;
    }

    public static int threeSumMulti(int[] arr, int target) throws InterruptedException {
        HashMap<Integer, Integer> frequencyMap = new HashMap<>();
        List<Integer> reverse = Arrays.stream(arr).boxed().collect(Collectors.toSet()).stream()
                .collect(Collectors.toList());

        // Collections.reverse(reverse);
        Collections.sort(reverse);

        if (arr.length == 3) {
            int sum = 0;
            for (int i = 0; i < arr.length; i++) {
                sum += arr[i];
            }

            return sum == target ? 1 : 0;
        }

        for (int i = 0; i < arr.length; i++) {
            frequencyMap.put(arr[i], frequencyMap.containsKey(arr[i]) ? frequencyMap.get(arr[i]) + 1 : 1);
        }

        int i = 0, j = 0;

        long multiplierCounter = 0L;

        System.out.println(frequencyMap.toString());

        while (i < reverse.size()) {
            // System.out.println(reverse.get(j) + " first");
            while (j < reverse.size()) {
                int pair = target - reverse.get(i) - reverse.get(j);
                System.out.println(pair + " " + reverse.get(i) + "<>" + reverse.get(j) + "<>");
                if (reverse.get(i) <= reverse.get(j) && reverse.get(j) <= pair && frequencyMap.containsKey(pair)
                        && pair + reverse.get(j) + reverse.get(i) == target) {
                    System.out.println(reverse.get(i) + "<a>" + reverse.get(j) + " " + pair);
                    /*
                     * Ahora comprueba los 3 pinches casos
                     * recuerda que inviertes el arreglo, so el j es el número mayor e i el menor
                     * we, comparas al reves
                     */

                    // **Este si fue hecho conb ayuda, la formula de los combinatronics no
                    // es un tema nativo de programación, repasarlo. Con brute force a cualquiera le
                    // queda,
                    // pero el tema de la calculación ahi si no */
                    int x = reverse.get(i);
                    int y = reverse.get(j);
                    int z = pair;
                    long fx = frequencyMap.get(x);
                    long fy = frequencyMap.get(y);
                    long fz = frequencyMap.get(pair);
                    if (x == y && y == z) {
                        multiplierCounter += (fx * (fx - 1) * (fx - 2) / 6) % 1_000_000_007;
                    } else if (x == y && y < z) {
                        multiplierCounter += ((fx * (fx - 1) / 2) * fz) % 1_000_000_007;
                    } else if (x < y && y == z) {
                        multiplierCounter += (fx * (fy * (fy - 1) / 2)) % 1_000_000_007;
                    } else {
                        multiplierCounter += (fx * fy * fz) % 1_000_000_007;
                    }

                    System.out.println("SUMA TOTAL " + multiplierCounter);
                    // multiplierCounter = 0;

                }

                j++;
            }
            i++;
            j = i;
        }

        return (int) multiplierCounter % 1_000_000_007;
    }

    /**
     * Aprendetelo para resolver los subsiguientes, sin letras son los más fáciles
     */
    public static List<List<Integer>> permute(int[] nums) {

        List<List<Integer>> result = new ArrayList<>();

        int[] recursiveCounterByLevel = new int[nums.length];
        Arrays.fill(recursiveCounterByLevel, 0);

        result.add(Arrays.stream(nums).boxed().collect(Collectors.toList()));

        /**
         * Algorithm:
         * 
         * The algorithm generates (n-1)! permutations of the first n-1 elements,
         * adjoining the last element to each of these.
         * This will generate all of the permutations that end with the last element.
         * If n is odd, swap the first and last element and if n is even,
         * then swap the ith element (i is the counter starting from 0)
         * and the last element and repeat the above algorithm till i is less than n.
         * In each iteration, the algorithm will produce
         * all the permutations that end with the current last element.
         */

        int i = 0;
        while (i < nums.length) {
            System.out.println("level " + i + " "
                    + Arrays.stream(recursiveCounterByLevel).boxed().collect(Collectors.toList()).toString());
            if (recursiveCounterByLevel[i] < i) {
                /** Comienza el swapping */
                if (i % 2 == 0) {
                    int elem = nums[i];
                    nums[i] = nums[0];
                    nums[0] = elem; /** Swapea el elemento 0 con el último */

                } else {
                    int elem = nums[recursiveCounterByLevel[i]];
                    nums[recursiveCounterByLevel[i]] = nums[i];
                    nums[i] = elem; /** Swapea el elemento 0 con el último */
                }

                result.add(Arrays.stream(nums).boxed().collect(Collectors.toList()));

                recursiveCounterByLevel[i]++;
                i = 0;

                // System.out.println(result.toString()+" asd");
            } else {
                recursiveCounterByLevel[i] = 0;
                i++;
                // System.out.println(result.toString()+" bsd");
            }
        }

        System.out.println(
                "level  " + Arrays.stream(recursiveCounterByLevel).boxed().collect(Collectors.toList()).toString());
        return result;

    }

    public static List<List<Integer>> subsetsWithDup(int[] nums) {
        Set<List<Integer>> output = new LinkedHashSet<>();
        int i = 0, j = 0;

        int[] recursiveC = new int[nums.length];
        Arrays.fill(recursiveC, 0);

        if (nums.length > 0)
            output.add(new ArrayList<>()); // se añade el vacio por regla

        while (i < nums.length) {
            List<Integer> pivotLister = new ArrayList<>();

            recursiveC[i] = 1;

            // pivotLister.add(nums[i]);

            System.out.println(Arrays.stream(recursiveC).boxed().collect(Collectors.toList()));
            while (j < nums.length) {
                // pivotLister.add(nums[i]);
                if (recursiveC[j] != 1) {
                    pivotLister.add(nums[j]);
                    // output.add(new ArrayList<>());
                    // output.add(new ArrayList<>(pivotLister));
                }
                // pivotLister.clear();
                j++;
            }

            output.add(new ArrayList<>(pivotLister));
            recursiveC[i] = 0;
            i++;
            j = 0;
        }

        return output.stream().collect(Collectors.toList());
    }

    public static String zigZagConversion(String s, int numRows) {
        if (numRows <= 1) {
            return s;
        }

        Queue<Character> helperPush = new LinkedList<>(
                s.chars().mapToObj(e -> Character.valueOf((char) e)).collect(Collectors.toList()));
        List<List<Character>> prueba = new ArrayList<>();

        String h = "";
        int i = 0, j = 0;
        String newStr = "";

        Set<Character> builderStr = new LinkedHashSet<>();
        List<Character> builderStrOdd = new ArrayList<>();

        int k = 0;

        int maxLength = Integer.MIN_VALUE;
        while (helperPush.size() > 1) {
            h += String.valueOf(helperPush.poll());

            if (h.length() == numRows) {
                List<Character> adder = h.chars().mapToObj(e -> (char) e).toList();
                maxLength = Math.max(maxLength, adder.size());
                prueba.add(adder);

                h = h.substring(h.length() - 1, h.length());
            }
        }

        h += String.valueOf(helperPush.poll());
        System.out.println(prueba.toString());

        List<Character> adder = h.chars().mapToObj(e -> (char) e).toList();
        prueba.add(adder);

        maxLength = Math.max(maxLength, adder.size());

        while (i < prueba.size()) {

            if (i % 2 != 0) {
                List<Character> piv = new ArrayList<>(prueba.get(i));
                // piv.set(0, '-');
                List<Character> piv2 = new ArrayList<>();

                if (piv.size() < maxLength) {
                    piv2.add('-');

                    // piv2.addAll(piv);
                }
                Collections.reverse(piv);

                piv2.addAll(piv);
                prueba.set(i, piv2);
            }

            i++;
        }

        System.out.println(prueba.toString());
        while (newStr.length() < s.length()) {
            List<Character> piv = new ArrayList<>(prueba.get(j));
            if (piv.size() > 0) {
                if (k % 2 != 0) {
                    builderStrOdd.add(piv.get(0));
                } else {
                    builderStr.add(piv.get(0));
                }
                // newStr+=piv.get(0);
                piv.remove(0); // remove the head
                prueba.set(j, piv);
            }

            j++;

            if (j >= prueba.size()) {
                newStr += k % 2 == 0 ? builderStr.stream().map(e -> String.valueOf(e)).collect(Collectors.joining())
                        : builderStrOdd.stream().map(e -> String.valueOf(e)).collect(Collectors.joining());
                // System.out.println(newStr);
                j = 0;
                k++;
                builderStr.clear();
                builderStrOdd.clear();
            }
        }

        return newStr.replace("-", "");
    }

    /* Del pdf llevo dia 1, 2, 4,8 */

    public static void printNumbers() {
        int[] counts = { 0, 0, 0, 0 };

        int i = counts.length - 1, j = counts.length - 1;

        while (counts[0] != 9) {
            System.out.println("i " + i);
            if (counts[i] == 9) {
                int aux = i;

                while (aux < counts.length) {
                    if (counts[aux] == 9) {
                        counts[aux] = 0;
                    }
                    aux++;
                }

                i--;
                counts[i]++;
            } else {
                counts[i]++;
                // System.out.println("YALXXX "+" "+counts[i]);
            }

            int auxPrint = 0;
            while (auxPrint < counts.length) {
                System.out.print(counts[auxPrint] + "<");
                auxPrint++;
            }

            System.out.println();
            i = counts.length - 1;
        }

    }

    public static int myAtoi(String s) {
        int result = 0;

        String concatter = "";
        Queue<Character> cleansedOutput = new LinkedList<>(
                s.trim().chars().mapToObj(r -> Character.valueOf((char) r)).collect(Collectors.toList()));
        int state = 0; // estado del automata

        Set<Character> charsProhibited = new HashSet<>(List.of('.', ' ', '+', '-')); // filtros para el estado dos

        Character current = null;

        while (cleansedOutput.size() > 0) {
            current = cleansedOutput.poll(); // obten la primera letra;

            if (Character.isLetter(current))
                break; // si detectas cualquier letra en el ciclo, ya no concatenes nada, lo cortas
            // if(charsProhibited.contains(current) && state >0) break;

            if (state == 0) { // en el primer estado tu esperas solo el -
                if (current == '-' || current == '+') {
                    concatter += String.valueOf(current); // concatena el simbolo de menos
                    state = 1; // como ya tienes un simbolo, pasa al siguiente estado
                } else if (Character.isDigit(current)) {
                    if (current != '0') {
                        concatter += String.valueOf(current); // si en el estado 0 encuentras un número, saltate a
                                                              // validar al estado 2
                        state = 2;
                    } else {
                        state = 1;
                    }
                } else {
                    concatter = "";
                    break; //
                }
            } else if (state == 1) { // si en el estado 1 validas puros ceros, omite esos leading ceros en el
                                     // cocatenado
                if (charsProhibited.contains(current)) {
                    concatter = "";
                    break;
                }
                if (current != '0' && Character.isDigit(current)) { // si es diferente a 0 y no es letra, saltate al 2
                                                                    // para validad números
                    concatter += String.valueOf(current);
                    state = 2;
                } else {
                    state = 1;
                }
            } else { // en este estado solo evaluas números
                if (charsProhibited.contains(current))
                    break;
                if (Character.isDigit(current))
                    concatter += String.valueOf(current);
            }

        }

        System.out.println(concatter + " resultado" + state + " " + current);
        // }

        if (concatter.length() == 0)
            return 0;
        if (state == 0 && Character.isLetter(current)) {
            concatter = "";
            return 0;
        }
        if (state == 1 && (charsProhibited.contains(current) || Character.isLetter(current)))
            return 0;

        int trueLength = concatter.charAt(0) == '-' || concatter.charAt(0) == '+'
                ? concatter.substring(1, concatter.length()).length()
                : concatter.length();

        if (trueLength > 10) {
            return Character.isDigit(concatter.charAt(0)) ? Integer.MAX_VALUE
                    : (concatter.charAt(0) == '-' ? Integer.MIN_VALUE : Integer.MAX_VALUE);
        }

        result = (int) Math.max(Integer.MIN_VALUE, Math.min(Integer.MAX_VALUE, Long.parseLong(concatter)));
        return result;
    }

    public static boolean isPalindrome(int x) {
        if (x < 0)
            return false;
        int limitCalculation = (int) Math.log10(x) + 1; // formula que te da las veces para dividir un número y sacar el
                                                        // decimal
        Double numberParsedForCalc = Double.valueOf((x));

        int i = 0;
        int position = 0;

        List<Integer> aux = new ArrayList<>();

        while (i < limitCalculation) {
            position = (int) Math.floor(numberParsedForCalc) % 10;
            numberParsedForCalc = numberParsedForCalc / 10;
            aux.add(position);
            i++;
        }

        /**
         * Que pendejez wey, si ya lo andas leyendo de atras jajaja, pero bueh, lei este
         * consejo de reddit, alch no lo pensé asi
         */
        List<Integer> reversedAux = new ArrayList<>(aux);
        Collections.reverse(reversedAux);

        return reversedAux.hashCode() == aux.hashCode() ? true : false;
    }

    public static int reverse(int x) {
        if (x == Integer.MAX_VALUE || x == Integer.MIN_VALUE)
            return 0;
        int limitCalculation = x > 0 ? (int) Math.log10(x) + 1 : (int) Math.log10(x * -1) + 1; // formula que te da las
                                                                                               // veces para dividir un
                                                                                               // número y sacar el
                                                                                               // decimal
        Double numberParsedForCalc = x > 0 ? Double.valueOf((x)) : Double.valueOf((x * -1));

        System.out.println(numberParsedForCalc);

        System.out.println(limitCalculation + " dddd");

        int i = 0;
        int position = 0;

        int concatter = 0;

        boolean limited = false;

        while (i < limitCalculation) {
            position = (int) Math.floor(numberParsedForCalc) % 10;

            numberParsedForCalc = numberParsedForCalc / 10;

            concatter = concatter * 10 + position;
            System.out.println(concatter);

            if (concatter > Integer.MAX_VALUE
                    || (concatter > (int) (Double.valueOf(Integer.MAX_VALUE) / 10 / 10 / 10 / 10 / 10) && i == 4
                            && limitCalculation == 10)) {
                // System.out.println("EHNTRA?");
                limited = true;
                break;
            }

            // aux.add(position);
            i++;
        }

        System.out.println(concatter);

        if (!limited && x < 0) {
            concatter = concatter * -1;
        }
        return limited ? 0 : concatter;
    }

    /** 0,0 , 0,1 - 0,2 */

    /**
     * De haber sabido que solo abarcaba las primeras letras, pero leetcode luego se
     * cotiza. Para ser un starterr que
     * lleva 2 meses de su vida leetcodeando, no está mal por más larga que alla
     * diseñado la solucion
     */
    public static String longestCommonPrefix(String[] strs) {
        HashMap<Integer, HashSet<String>> prefixCounter = new HashMap<>();
        HashMap<String, Integer> charAcrossCounter = new HashMap<>();
        boolean sequential = false;
        String helper = "";

        Optional<Integer> minSize = Arrays.asList(strs).stream().map(e -> e.length()).min(Comparator.naturalOrder());

        if (!minSize.isPresent())
            return ""; // casos vacios

        int subsSlice = minSize.get();

        for (String stringToSplit : strs) {
            /// substrings.add();
            String checkerPerElem = stringToSplit.substring(0, subsSlice);
            int temp = 0;

            while (temp < checkerPerElem.length()) {
                // System.out.println(checkerPerElem.charAt(temp));
                if (!prefixCounter.containsKey(temp)) {
                    HashSet<String> tempArr = new HashSet<>();
                    tempArr.add(String.valueOf(checkerPerElem.charAt(temp)));
                    prefixCounter.put(temp, new HashSet<>(tempArr));
                } else {
                    // System.out.println(checkerPerElem.charAt(temp)+" existe "+temp);
                    HashSet<String> tempArr = prefixCounter.get(temp);
                    tempArr.add(String.valueOf(checkerPerElem.charAt(temp)));
                    prefixCounter.put(temp, new HashSet<>(tempArr));
                }

                temp++;
            }

        }

        System.out.println(prefixCounter.toString());

        Map<Integer, String> maximumLength = new TreeMap<>(Comparator.reverseOrder());

        int c = 0;

        while (c < prefixCounter.size()) {
            if (prefixCounter.get(c).size() == 1) {
                if (!sequential && c == 0) {
                    sequential = true;
                }
                if (sequential) {
                    helper += prefixCounter.get(c).iterator().next(); // como es set para evitar duplicados usate el
                                                                      // next del iterator
                }

            } else {
                if (sequential) {
                    maximumLength.putIfAbsent(helper.length(), helper);
                    helper = "";
                    sequential = false;
                    break; // porque solo debe corresponder a las primerasd letras péro leetcode siempre
                           // pone sus testcases ambiguos
                }
            }
            c++;
        }

        /**
         * DCon este snippet comentado, System.out.println(longestCommonPrefix(new
         * String[]{"flower","fkow"})+" respuesta"); este test case
         * te devuelve el true, lo pensé para indexes en comun
         */
        /*
         * for (int i = 0; i < prefixCounter.size(); i++) {
         * if( prefixCounter.get(i).size() ==1){
         * if(!sequential){
         * sequential = true;
         * }
         * if(sequential){
         * helper += prefixCounter.get(i).iterator().next(); // como es set para evitar
         * duplicados usate el next del iterator
         * }
         * 
         * }else{
         * if(sequential){
         * maximumLength.putIfAbsent(helper.length(), helper);
         * helper="";
         * sequential = false;
         * }
         * }
         * }
         */
        maximumLength.putIfAbsent(helper.length(), helper);
        prefixCounter.clear();

        System.out.println(maximumLength.toString());
        return maximumLength.size() > 0 ? (maximumLength.entrySet().iterator().next().getValue().isEmpty() ? ""
                : maximumLength.entrySet().iterator().next().getValue()) : "";

    }

    public static int[] twoSum(int[] nums, int target) {
        int[] results = new int[2];

        HashMap<Integer, Queue<Integer>> indexSaver = new HashMap<>();
        Set<Integer> complementSearcher = new HashSet<>(Arrays.stream(nums).boxed().collect(Collectors.toSet()));

        for (int i = 0; i < nums.length; i++) {
            Queue<Integer> tempArr = new LinkedList<>();
            if (!indexSaver.containsKey(nums[i])) {
                tempArr.add(i);
                indexSaver.put(nums[i], new LinkedList<>(tempArr));
            } else {
                tempArr.addAll(indexSaver.get(nums[i]));
                tempArr.add(i);
                indexSaver.put(nums[i], new LinkedList<>(tempArr));
            }
        }

        // System.out.println(indexSaver.toString());

        int i = 0;
        while (i < nums.length) {
            int complement = target - nums[i];

            if (complement != nums[i]) {
                if (complementSearcher.contains(complement)) {
                    results[0] = i;
                    results[1] = indexSaver.get(complement).poll(); // remueve uno o más indices dependiendo de si hay
                                                                    // números duplicados
                    break;
                }
            } else {
                if (complementSearcher.contains(complement) && indexSaver.get(complement).size() > 1) {
                    results[0] = indexSaver.get(complement).poll(); // remueve uno o más indices dependiendo de si hay
                                                                    // números duplicados
                    results[1] = indexSaver.get(complement).poll();
                    break;
                }
            }

            i++;
        }

        /*
         * for (Integer integer : results) {
         * System.out.print(integer);
         * }
         * System.out.println();
         */
        // System.out.println(complementSearcher.toString());

        return results;
    }

    static class ListNode {
        int val;
        ListNode next;

        ListNode() {
        }

        ListNode(int val) {
            this.val = val;
        }

        ListNode(int val, ListNode next) {
            this.val = val;
            this.next = next;
        }

        public String toString(){
            return "OBJ "+val;
        }

    }

    public static ListNode addTwoNumbers(ListNode l1, ListNode l2) {
        ListNode resultSum = null;
        ListNode c1 = l1,c2 = l2;
        BigInteger numOne = new BigInteger("0"),numTwo = new BigInteger("0");
        Stack<Integer> c1Reverse = new Stack(),c2Reverse = new Stack();

        while (c1 !=null) {
            c1Reverse.push(c1.val);
            c1 = c1.next;
        }

        while(c2 !=null){
            c2Reverse.push(c2.val);
            c2 = c2.next;
        }

    //    System.out.println("ESTADO DE LAS PILAS "+c1Reverse.toString()+" "+c2Reverse.toString());

        while(c1Reverse.size() >0){
          numOne = new BigInteger(String.valueOf(numOne).concat(String.valueOf(c1Reverse.pop())));
        }

        while(c2Reverse.size() >0){
          numTwo = new BigInteger(String.valueOf(numTwo).concat(String.valueOf(c2Reverse.pop())) );
         
        }

        String resultToTransform = String.valueOf(numOne.add(numTwo) ); // de 807 tienes que convertirlo a 708, como lo logras

        int cc = 0;

        
        while(cc<resultToTransform.length()){
            ListNode tempNode = new ListNode();
            int charNum = resultToTransform.charAt(cc) - '0'; //conviertelo a número
            tempNode.val = charNum;

            if(resultSum ==null){
                tempNode.next = null;
                resultSum = tempNode;
            }else{
                tempNode.next = resultSum;
                resultSum = tempNode;
            }
            cc++;
        }

        return resultSum;
    }

    public static double findMedianSortedArrays(int[] nums1, int[] nums2) {
        double result = 0d;

        int[] generalNums = new int[nums1.length+nums2.length];

        /**Popula el arreglo */
        for (int i = 0; i < nums1.length; i++) {
            generalNums[i] = nums1[i];
        }

        int iterFollow = nums1.length;
        for (int i = 0; i < nums2.length; i++) {
            generalNums[iterFollow] = nums2[i];
            iterFollow++;
        }


        /**Sortea el arreglo en orden natural */
        Arrays.sort(generalNums);

        int middleOne = 0,middleTwo = 0;

      //  System.out.println(generalNums.length);

        if(generalNums.length %2 == 0){
            middleOne = (generalNums.length-1) /2;
            middleTwo = middleOne + 1; //obtenemos el número de la izquierda
            System.out.println(generalNums[middleOne] +"<<"+ generalNums[middleTwo]);
            result = (double) (((double) generalNums[middleOne] + (double)generalNums[middleTwo]) / 2);
        }else{
            middleOne = (generalNums.length-1) /2;
            result = generalNums[middleOne];
        }

        
        return result;
    }

    public static String convert(String s, int numRows) throws InterruptedException {
        String finalStr = "";

        char[][] lettersZigZag = null;// n+1 * n para la fila extra pa

        System.out.println();

        lettersZigZag = new char[s.length()][numRows];
      
        if (numRows <= 1 || s.length()<=numRows) {
            return s;
        }
        
        int globalChar = 0;


        /**La idea con este approach es concatenar esto =
         * 
         * 
         * PINALSIGYAHRPI
            [
                [P, A, Y, P], 
                [I, L, A, -], 
                [-, S, H, I], 
                [N, I, R, -], 
                [-, G, -, -]
            ]
         */
      System.out.println("tamaño del string "+s.length());
        for (int k = 0; k < lettersZigZag.length; k++) {
        //    System.out.println(k+" fase del cponmcatenao");
            globalChar=k>0 ? globalChar-=1:0; // para poner la letra del último ciclo en la proxima fila invertidamente
            if(k %2 ==0){
                for (int k2 = 0; k2 < lettersZigZag[0].length; k2++) {
                    if(globalChar<s.length()){
                        lettersZigZag[k][k2] = k2 == 0 && k>0 ? '-' : s.charAt(globalChar);
                        globalChar++;
                    }else{
                        lettersZigZag[k][k2] = '-';
                    }
                }
            }else{
          //      System.out.println(k+" fase del cponmcatenao");
                for (int k3 = lettersZigZag[0].length-1; k3 >=0; k3--) {
                    if(globalChar<s.length()){
                        lettersZigZag[k][k3] =  k3 == lettersZigZag[0].length-1  && k>0? '-' : s.charAt(globalChar);
                        globalChar++;
                    }else{
                        lettersZigZag[k][k3] = '-';
                    }
                }
            }
        }

        System.out.println(globalChar+" "+lettersZigZag.length+" "+lettersZigZag[0].length);

        for (int k = 0; k < lettersZigZag.length; k++) {
            for (int k2 = 0; k2 < lettersZigZag[0].length; k2++) {
                System.out.print(lettersZigZag[k][k2]);
            }
            System.out.println();
        }  

      //  System.out.println();System.out.println();
        for (int k = 0; k < lettersZigZag[0].length; k++) {
           for (int k2 = 0; k2 < lettersZigZag.length; k2++) {
                if(lettersZigZag[k2][k] != '-'){
                //    System.out.print(lettersZigZag[k2][k]);
                    finalStr+=String.valueOf(lettersZigZag[k2][k]);
                }
           }
        }

        return finalStr;

    }

    public static int longestPalindromeV2(String[] words) throws InterruptedException {
        int i = 0;
        StringBuilder c1 = new StringBuilder();
        StringBuilder c2 = new StringBuilder();

        Map<String,Integer> concurrencyCounter = new LinkedHashMap<>();

        for (int j = 0; j < words.length; j++) {
            if(!concurrencyCounter.containsKey(words[j])){
                concurrencyCounter.put(words[j], 1);
            }else{
                concurrencyCounter.put(words[j],concurrencyCounter.get(words[j])+1);
            }
        }



        System.out.println(concurrencyCounter.toString());

        Iterator<Entry<String,Integer>> iter = concurrencyCounter.entrySet().iterator();

        Map<Integer,String> maxLength = new TreeMap<>(Comparator.reverseOrder());

        int impairLast = 0,originalQ = 0;
        String lastRepCh = "";
        while(iter.hasNext()){
            Entry<String,Integer> subsChecker = iter.next();
            StringBuilder reverseHelper = new StringBuilder();
            String reverse = reverseHelper.append(subsChecker.getKey()).reverse().toString(); //reversa el char

            if(subsChecker.getValue() >0){ // para que en los ciclos no tome pares de más
                if(subsChecker.getValue() ==1 && concurrencyCounter.containsKey(reverse) && reverse.equals(subsChecker.getKey())){ // si solo corresponde a uno
                    c1.insert(c1.length(), reverse);
                    lastRepCh = reverse;
                    impairLast = 0;
                   // originalQ = subsChecker.getValue();
                    concurrencyCounter.put(reverse, 0);
                    concurrencyCounter.put(subsChecker.getKey(), 0);
                    
                    maxLength.put(c1.length()+c2.length(), c1.toString().concat(c2.toString()).toString());
                
                    /*Reinicia tus concatters, porque cuando encuentras un par que solo tiene uno, no puedes añadir más repetidos de par en medio */
                    c1 = new StringBuilder();
                    c2= new StringBuilder();
                    
                }else if(subsChecker.getValue() ==1 && concurrencyCounter.containsKey(reverse) && !reverse.equals(subsChecker.getKey())){ /**Pares distintos */
                    c1.insert(0, subsChecker.getKey());
                    c2.insert(c2.length(), reverse);

                    concurrencyCounter.put(reverse, 0);
                    concurrencyCounter.put(subsChecker.getKey(), 0);

                    maxLength.put(c1.length()+c2.length(), c1.toString().concat(c2.toString()).toString()); // nota mental del builder, no uses append porque afecta al objeto, siempre convierte el to string 
                    
                }else if(subsChecker.getValue() > 1 && concurrencyCounter.containsKey(reverse) && !reverse.equals(subsChecker.getKey())){
                    /**Puedes tener 3 lc y 10 cl's, en base al menor valor, iras convirtiendo al cero y despues del ciclo, setear en 0 el cl
                     * puesto que el residuo ya no lo usarás, nisiquiera en medio
                     */
                    int minIter = Math.min(concurrencyCounter.get(reverse), subsChecker.getValue()); //3

     
                    c1.insert(0, subsChecker.getKey().repeat(minIter));
                    c2.insert(c2.length(), reverse.repeat(minIter));

                    concurrencyCounter.put(reverse, 0);
                    concurrencyCounter.put(subsChecker.getKey(), 0);

                    maxLength.put(c1.length()+c2.length(), c1.toString().concat(c2.toString()).toString());

                }else if(subsChecker.getValue() > 1 && concurrencyCounter.containsKey(reverse) && reverse.equals(subsChecker.getKey())){
                    /**En este escenario, hay pares repetidos, teniendo que calcular cuantos al lado y el residuo central, pero también
                     * sobre eso, no repetir
                     */

                    int maxPairRepeated = (int)subsChecker.getValue() / 2;
                    int res = (int)subsChecker.getValue() % 2;
                    lastRepCh = reverse;
                    impairLast = res;

                    c1.insert(c1.length(), subsChecker.getKey().repeat(maxPairRepeated));
                    c2.insert(0, reverse.repeat(maxPairRepeated));
                    
                    concurrencyCounter.put(reverse, 0);
                    concurrencyCounter.put(subsChecker.getKey(), 0);

                    maxLength.put(c1.length()+c2.length(), c1.toString().concat(c2.toString()).toString());
                }

            }
            
            
        }
        System.out.println(maxLength.toString()+" generados "+impairLast);
        

        if(maxLength.size() ==0) return 0;
        if(impairLast == 0) return maxLength.entrySet().iterator().next().getKey();
        Entry<Integer,String> pair = maxLength.entrySet().iterator().next();
        String finalConcat = pair.getValue().substring(0, (pair.getValue().length())/2);
        return String.valueOf(finalConcat+lastRepCh.repeat(impairLast)+pair.getValue().substring((pair.getValue().length())/2,pair.getValue().length()) ).length();
    }


    public static int numberOfSubarrays(int[] nums, int k) {
        int one = 0, two = 1; // pointers to move

     //   Queue<Integer> odds = new LinkedList<>();

        int aux = 0;
        int counterK = 0,counterTotal=0,counterPerCycle = 0;


        while (one <= nums.length-2 ) {
            System.out.println("ciclos");
            if(nums[one] %2 !=0){
                //oddInterceptors.add(one);
                counterK++;
            }
            //counterPerCycle++;

            while(two<nums.length){

                if(counterK ==k){
                    counterTotal++;
                }else if(counterK>k){
                    counterTotal+=counterPerCycle;
                    counterTotal++;
                    break;
                }

                if(nums[two] %2 !=0){
                  //  oddInterceptors.add(two);
                    counterK++;
                }

                counterPerCycle++;

                System.out.println(nums[two]+" "+counterPerCycle);

                two++;
            }

          /*   if(counterK ==k){
                counterTotal++;
            } */
            counterPerCycle=0;

            
            one++;
            two = one +1;
             
        }


        return counterTotal;
    }


    public static int fourSumCount(int[] nums1, int[] nums2, int[] nums3, int[] nums4) {
    
        int toupleCount=0;

      /*   Arrays.sort(nums1);
        Arrays.sort(nums2);
        Arrays.sort(nums3);
        Arrays.sort(nums4); */

        HashMap<Integer,Integer> lastNum = new HashMap<>();
        HashMap<Integer,Integer> preNum = new HashMap<>();

        for (int i = 0; i < nums4.length; i++) {
            lastNum.put(nums4[i], lastNum.containsKey(nums4[i]) ?lastNum.get(nums4[i])+1 : 1 );
        }

        for (int i = 0; i < nums3.length; i++) {
            preNum.put(nums3[i], preNum.containsKey(nums3[i]) ?preNum.get(nums3[i])+1 : 1 );
        }

        System.out.println(lastNum.toString());
      //  System.out.println(preNum.toString());

        for (int i = 0; i < nums1.length; i++) {
           // System.out.println(nums1[i]+" i");
            for (int j = 0; j < nums2.length; j++) {
              //  System.out.println(nums2[j]+" j");
              /*   int preSum = nums1[i] + nums2[j];
                if(preNum.containsKey(preSum*-1) ){
                    toupleCount+=preNum.get(preSum*-1);
                    preSum +=(preSum*-1);
                    if(lastNum.containsKey(preSum * -1)){
                        toupleCount+=lastNum.get(preSum*-1);
                    }
                } */
                 for (int k = 0; k < nums3.length; k++) {
                    int preSum = nums1[i] + nums2[j] + nums3[k];
                    if(lastNum.containsKey(preSum * -1)){
                        toupleCount+=lastNum.get(preSum*-1);
                    }
                } 
            }
        }

        return toupleCount;
    }

    public static void spiralOrder(int[][] matrix) throws InterruptedException{
        int top = 0,right = matrix[0].length,bottom = matrix.length-1,left=0;

        int i=0,j=0;

        while(top <= bottom && left<=right){
            while(j<top){
                System.out.println(matrix[i][j]);
                j++;
            }
            top++;
           // j=top;
            
            while(i<right){
                System.out.println(matrix[i][j]);
                i++;
            } right--;

          //  i=right;

            while(j>=bottom){
                System.out.println(matrix[i][j]);
                j--;
            }

           // j=bottom;

           // Print bottom row from right to left (if exists)
            if (top <= bottom) {
                 while(j>=left){
                    System.out.println(matrix[i][j]);
                    if(i==left) break;
                    j--;
                }
                bottom--;
            }

            if(left<=right ){
                System.out.println("SII "+left+" "+i);
                while(i>=left){
                    System.out.println(matrix[i][j]);
                    if(i==left) break;
                    i--;
                }
                left++;
            }
            
        }
    }


    public static int uniquePathsWithObstacles(int[][] obstacleGrid) {
        
        Stack<HashMap<Integer,Integer>> visiting = new Stack<>();

        Set<HashMap<Integer,Integer>> visited = new HashSet<>();

        int i=0,j=0; int genCounter = 0;

        /**Arriba, izquierda, derecha, abajo */
      /*   int[] row = {-1,0,0,1};
        int[] col = {0,-1,1,0}; */

        int[] row = { -1, 0, 0, 1 };
        int[] col = { 0, -1, 1, 0 };

        if(obstacleGrid[0][0] == 1 || obstacleGrid[obstacleGrid.length-1][obstacleGrid[0].length-1] == 1) return 0; // obstaculo en la entrada lol
       // if(obstacleGrid.length ==1 && (obstacleGrid[0][0] ==1 || obstacleGrid[obstacleGrid.length-1][obstacleGrid[0].length-1] == 1)) return0;

        if(obstacleGrid.length ==1){
            Set<Integer> checker = new HashSet<>();
            checker.addAll(Arrays.stream(obstacleGrid[0]).mapToObj(e->(int)e).collect(Collectors.toSet()));

            return checker.size()==1 ? 1:0;
        }

        while(i<obstacleGrid.length){
            while (j<obstacleGrid[0].length) {
                HashMap<Integer,Integer> coord = new HashMap<>();
                coord.put(i, j);
                if(obstacleGrid[i][j] != 1 && !visited.contains(coord)){
                    visiting.push(new HashMap<>(coord)); //coordenada actual como visitada
               //     visiting.push(new HashMap<>(coord));
                    
                    while (visiting.size()>0) {
                        HashMap<Integer,Integer> inst = visiting.pop();
                        Entry<Integer,Integer> opc =  inst.entrySet().iterator().next();
                        //HashMap<Integer,Integer> 
                       
                        
                        if(opc.getKey() == obstacleGrid.length-1 && opc.getValue() == obstacleGrid[0].length-1){
                            System.out.println("YEEEEEEEEEEEEEEEEEE");
                            visiting.clear();
                            visited.clear();
                            genCounter++;
                        }else{
                            visited.add(inst);
                            if(obstacleGrid[opc.getKey()][opc.getValue()] !=1){
                                for (int k = 0; k < col.length; k++) {
                                    int rowCalc = opc.getKey()+row[k];
                                    int colCalc = opc.getValue()+col[k];

                                    System.out.println(rowCalc+"-"+colCalc);

                                    if(rowCalc>=0 && rowCalc<obstacleGrid.length && colCalc>=0 && colCalc<obstacleGrid[0].length){
                                        if(obstacleGrid[rowCalc][colCalc] ==0){
                                            HashMap<Integer,Integer> crCoord = new HashMap<>();
                                            crCoord.put(rowCalc, colCalc);
                                            if(!visited.contains(crCoord)){
                                                visiting.add(new HashMap<>(crCoord));
                                            }
                                        }
                                        
                                    }
                                }
                            }
                        }
                        System.out.println(visiting.toString());
                    }
                     visited.add(coord);
                }
                j++;
            }

        
            i++;
        }

        System.out.println(genCounter+" asDASDASDSADSA");
        return genCounter;
    }

    public static int maxProfit(int[] prices) {
        TreeMap<Integer,TreeSet<Integer>> bestProf = new TreeMap<>(Comparator.naturalOrder());
        for (int j2 = 0; j2 < prices.length; j2++) {
            int mainC = prices[j2];
            for (int k = j2+1; k < prices.length; k++) {
             //   System.out.println(prices[k]+" d");
                if(mainC<prices[k]){
                //    System.out.println("NEXT "+(Math.abs((mainC - prices[k])  *-1))+" >>"+mainC);
                    if(!bestProf.containsKey(mainC)){
                        TreeSet<Integer> aux = new TreeSet<>(Comparator.reverseOrder());
                        aux.add(Math.max(Math.abs((mainC - prices[k])), Math.abs((mainC - prices[k])) * -1));
                        bestProf.put(mainC, new TreeSet<>(aux));
                    }else{
                        TreeSet<Integer> aux = bestProf.get(mainC);
                        aux.add(Math.max(Math.abs((mainC - prices[k])), Math.abs((mainC - prices[k])) * -1));
                        bestProf.put(mainC, new TreeSet<>(aux));
                    }
                 //   System.out.println(bestProf.toString());
                }else{
                  //  System.out.println("NO");
                    if(bestProf.containsKey(mainC)){
                        TreeSet<Integer> aux = bestProf.get(mainC);
                        aux.add(0);
                        bestProf.put(mainC, new TreeSet<>(aux));
                    }else{
                        TreeSet<Integer> aux = new TreeSet<>(Comparator.reverseOrder());
                        aux.add(0);
                        bestProf.put(mainC, new TreeSet<>(aux));
                    }
                }
                
               
               // System.out.println();
            }
        }

     //   System.out.println(bestProf.toString());

        return bestProf.size()> 0 ? bestProf.entrySet().iterator().next().getValue().iterator().next():0;
    }

    public static void main(String[] args) throws InterruptedException {

         System.out.println(maxProfit(new int[]{7,1,5,3,6,4}));

         System.out.println(maxProfit(new int[]{7,6,4,3,1})); 
         System.out.println(maxProfit(new int[]{2,4,1}));
/*         System.out.println(uniquePathsWithObstacles(new int[][]{{0,0,0},{0,1,0},{0,0,0}})+" otro bfs we");

        System.out.println(uniquePathsWithObstacles(new int[][]{{0,1},{0,0}})+" otro bfs we");

         System.out.println(uniquePathsWithObstacles(new int[][]{{0,0},{1,0}})+" otro bfs we"); */

       //  System.out.println(uniquePathsWithObstacles(new int[][]{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,1,0},{0,0,0,0}})+" otro bfs we");
      //  spiralOrder(new int[][]{{1,2},{4,5}});
    //    spiralOrder(new int[][]{{1,2,3},{4,5,6},{7,8,9}});

     //   spiralOrder(new int[][]{{1,2,3,4},{4,5,6,7},{7,8,9,10}});
 /* System.out.println(fourSumCount(new int[]{1,2} ,
            new int[]{-2,-1},
            new int[]{-1,2},
            new int[]{0,2}));
            
        System.out.println(fourSumCount(new int[]{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2}, 
            new int[]{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2}, 
            new int[]{-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2},
            new int[]{-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2}));
 */
     /*    System.out.println(fourSumCount(new int[]{59,-2,34,81,-91,4,78,60,96,-25,-19,33,-39,23,-60,4,-65,-9,77,-52,54,5,-64,8,73,-52,49,-6,-54,34,60,99,-4,-54,62,63,-7,-1,-60,43,-46,-61,14,-86,-23,83,-92,40,-32,-27,-70,81,6,52,64,76,89,48,21,19,-61,38,-88,36,75,-83,44,-33,-81,-50,-7,45,79,80,-1,74,40,-80,-61,46,39,-50,-27,-56,35,97,31,23,90,-30,-84,-61,-61,22,40,-69,-38,0,-65,-92,36,0,-7,-32,-36,21,-39,-74,93,9,-33,22,-81,-36,49,81,43,-84,71,40,78,50,-1,-48,-98,99,-70,36,3,-13,77,17,-32,4,39,-7,-25,-90,-7,15,72,52,26,100,-83,7,0,13,-41,-2,-34,-87,40,-75,2,36,44,74,-68,-96,-63,16,-89,29,12,-24,-74,-50,-70,8,34,5,61,-11,58,65,29,-19,-61,-14,-51,-92,35,-59,-75,-38,51,-66,41,-93,-88,-93,-85,-2},
            new int[]{98,54,-33,30,92,13,100,-34,34,24,-12,86,87,-78,31,22,23,-39,-18,41,87,-4,42,41,-43,16,-13,-39,26,99,15,72,28,93,95,12,48,-5,74,-78,57,-39,-29,-42,54,-23,35,43,-17,-58,-32,-67,15,-51,62,-79,49,1,-18,-45,71,61,-96,-78,46,-10,-76,69,-98,-93,-69,26,22,19,65,96,-87,-44,-67,-47,88,-51,-79,-81,22,97,62,7,-83,-57,-53,51,38,27,27,73,-9,1,-33,56,87,-3,28,-47,-83,-88,-56,-46,78,-7,-69,-18,-74,-20,-43,-21,13,9,81,-88,-27,-16,89,-30,3,41,-76,-45,64,-26,-26,-58,-86,-34,56,-97,75,21,-55,-63,-60,-10,-61,55,-3,-20,-99,-83,-54,-2,-57,-16,12,-92,-42,-70,97,-78,20,69,10,24,-90,-52,-65,98,15,-73,-25,-48,-98,0,86,59,40,67,-41,-72,95,-11,50,-41,77,16,54,3,49,-35,-90,12,17,80,100,-31},
            new int[]{59,87,-37,15,-33,-53,22,7,6,-12,65,90,-10,-65,40,-7,20,-76,5,-76,32,67,-82,-26,15,-34,-55,65,-79,98,-99,85,-32,22,43,13,-25,32,43,30,8,23,-21,82,98,-61,2,77,-24,12,60,39,-51,1,-95,-38,49,-1,-37,92,-30,90,24,-75,45,-11,-66,-11,-62,52,64,-17,-100,-81,-41,-28,80,-19,-78,-85,-99,55,33,72,7,-23,64,7,-26,2,-46,42,-54,15,-84,41,27,-87,-6,-82,35,100,-57,-99,-72,65,31,-53,13,-31,-46,65,-17,-88,34,68,-44,-48,-100,21,-82,70,-28,13,-47,87,15,-31,-72,-59,-69,-94,87,-30,84,-63,-11,-32,76,-59,29,-82,-100,54,49,-18,29,31,-65,70,-51,-56,7,59,-69,-35,82,40,79,-62,-96,29,-54,-69,-48,58,-4,73,-73,31,42,-1,-56,-7,-24,-63,83,-48,84,12,22,-20,16,94,-46,-44,-63,40,14,-35,-16,14,100,-28},
            new int[]{-91,38,-11,-89,95,-66,-87,85,25,50,-18,76,20,-37,-3,-67,-76,85,73,-3,-6,14,-60,-85,75,-46,6,-46,88,-78,0,-90,36,32,76,-7,86,-57,-30,-34,-56,-99,60,-9,60,-19,-21,-73,44,61,23,75,61,45,60,59,69,47,59,54,55,-79,-92,12,24,-62,75,-88,23,39,19,-74,43,-75,-20,18,-56,52,2,-80,-59,-9,44,42,-7,-6,12,-3,-12,85,-66,0,0,-80,50,88,-79,-77,43,22,59,-5,97,-28,-2,80,47,-41,53,-39,-41,30,-10,76,-3,32,7,99,12,79,8,-2,87,44,73,-46,89,1,-47,92,34,53,33,42,-95,93,-3,77,-14,64,-27,10,80,40,-77,51,-15,-2,34,-54,-19,-11,-41,65,-36,97,53,-60,93,-81,58,-27,-35,-97,67,26,-5,91,6,-62,-50,-9,-22,-17,-89,-37,-14,-22,10,77,-67,3,10,32,14,-94,-34,33,49,42,-4,-52,-22,0}));
 */
/*         System.out.println(numberOfSubarrays(new int[]{2,2,2,1,2,2,1,2,2,2}, 2)+" RESPUESTA");
        System.out.println(numberOfSubarrays(new int[]{1,1,2,1,1}, 3)+" RESPUESTA"); */
        /* System.out.println(numberOfSubarrays(new int[]{1,1,2,1,1}, 3)+" RESPUESTA");

        System.out.println(numberOfSubarrays(new int[]{2,4,6}, 1)+" RESPUESTA"); */

        
/* 
        System.out.println(longestPalindromeV2(new String[]{"lc","cl","gg"}));

        System.out.println(longestPalindromeV2(new String[]{"ab","ty","yt","lc","cl","ab"}));

        System.out.println(longestPalindromeV2(new String[]{"cc","ll","xx"}));

        System.out.println(longestPalindromeV2(new String[]{"lc","cl","gg","gg","gg","gg","gg","lc","lc","cl","cl","cl","cl","cl","cl","cl","cl","cl","cl"})); // PARA PROBAR TERCER CASO DE USO


        System.out.println("LOS PRIMEROS WE");
        System.out.println(longestPalindromeV2(new String[]{"dd","aa","bb","dd","aa","dd","bb","dd","aa","cc","bb","cc","dd","cc"}));
        System.out.println(longestPalindromeV2(new String[]{"lc","cl","gg","ak","kk"}));
        System.out.println(longestPalindromeV2(new String[]{"ab","ty","yt","lc","cl","ab","ba"}));
        System.out.println(longestPalindromeV2(new String[]{"qo","fo","fq","qf","fo","ff","qq","qf","of","of","oo","of","of","qf","qf","of"}));
        System.out.println(longestPalindromeV2(new String[]{"io","io"}));

        System.out.println(longestPalindromeV2(new String[]{"ga","ac","aa","ag","gc","cg","aa","ac","cg","ga","ga","gg","cg","ca","cg","gg","ca","ag","cc","ag","aa","cg","gg"}));
        System.out.println(longestPalindromeV2(new String[]{"ll","lb","bb","bx","xx","lx","xx","lx","ll","xb","bx","lb","bb","lb","bl","bb","bx","xl","lb","xx"}));
        System.out.println(longestPalindromeV2(new String[]{"qw","rr","ll","vv","iw","wq","cc","wi","jj","iw","pp","iw","mm","ss","bb","oo","wi","dd","wq","ff","qi","qw","qi","qi","zz","wq","iw","wi","qq","qw","wi","hh","qi","pp","vv","wi","wq","wi","wi","wi","iw","qi","bb","qw","qi","rr"}));
 */
        /*         
        System.out.println(convert("PAYPALISHIRING",4)+" RERSULTADO");
        System.out.println(convert("PAYPALISHIRING",3)+" RERSULTADO");
        System.out.println(convert("ABCDE",2)+" RERSULTADO");

        System.out.println(convert("ABCDEF",2)+" RERSULTADO");

        System.out.println(convert("PAYPALISHIRING",2)+" RERSULTADO"); */
        //System.out.println(findMedianSortedArrays(new int[]{1,2}, new int[]{3,4}));
/* 
        ListNode c1 = new ListNode();

        c1.val = 2;
        c1.next = new ListNode();
        c1.next.val = 4;
        c1.next.next = new ListNode();
        c1.next.next.val = 3; */

     /*    ListNode c2 = c1;
        while (c2 != null) {
            System.out.println(c2.toString());
            // System.out.println(c2.next.toString());
            c2 = c2.next;
        }
 */
  /*       ListNode c3 = new ListNode();

        c3.val = 5;
        c3.next = new ListNode();
        c3.next.val = 6;
        c3.next.next = new ListNode();
        c3.next.next.val = 4; */

       /*  ListNode c4 = c3;
        while (c4 != null) {
            System.out.println(c3.toString());
            c4 = c4.next;
        }
 */

/*        ListNode temp = addTwoNumbers(c1, c3);

       while (temp !=null) {
        System.out.println(temp.val);
        temp = temp.next;
       }

       BigInteger r1 = new BigInteger("1000000000000000000000000000001");
       BigInteger r2 = new BigInteger("1000000000000000000000000000001");

      BigInteger res1 = r1.add(r2);

       System.out.println(r1+" "+r2+" "+r1.add(r2)); */

        /** AQUI ESTÁ LÑA APLICACION DE ESTOS METODOS PARA PODERFT AOCIKASOFO O */

        /*
         * System.out.println(twoSum(new int[]{2,7,11,15}, 9));
         * System.out.println(twoSum(new int[]{3,2,4}, 6));
         * System.out.println(twoSum(new int[]{3,3}, 6));
         */
        // System.out.println(longestCommonPrefix(new
        // String[]{"flower","flow","flight"}));

        /*
         * System.out.println(longestCommonPrefix(new
         * String[]{"flower","flow","flight"})+" respuesta");
         * System.out.println(longestCommonPrefix(new
         * String[]{"rowdynas","rowpynas","rowllbynas"})+" respuesta");
         * System.out.println(longestCommonPrefix(new String[]{"a", "aca", "accb",
         * "b"})+" respuesta");
         * System.out.println(longestCommonPrefix(new
         * String[]{"iflower","flow","flight"})+" RESPUESTA");
         * System.out.println(longestCommonPrefix(new String[]{"c","acc","ccc"}));
         * System.out.println(longestCommonPrefix(new
         * String[]{"reflower","flow","flight"})+" respuesta");
         * 
         * System.out.println(longestCommonPrefix(new
         * String[]{"dog","racecar","car"})+" respuesta");
         */
        /*
         * System.out.println(longestCommonPrefix(new
         * String[]{"flower","flow","flight"}));
         * System.out.println(longestCommonPrefix(new String[]{"a"})+" respuesta");
         * System.out.println(longestCommonPrefix(new
         * String[]{"flower","fkow"})+" respuesta");
         * System.out.println(longestCommonPrefix(new
         * String[]{"aca","cba"})+" respuesta");
         */
        // System.out.println(longestCommonPrefix(new
        // String[]{"radosadasdadasdasdasdasdasdadasdasdasaddddddddddddddg","racasadasdasdasdasdasdasdasdcxcvzxccsddssczxcsdfcsdcsdcsdcsdcddcdscecar","racar","raracacar"}));
        /*
         * System.out.println((int)(Double.valueOf(Integer.MAX_VALUE)/10 /10/10/10/10
         * )+" "+Integer.MIN_VALUE+" >>>"+Integer.MAX_VALUE);
         * System.out.println(reverse(1563847412));
         */
        /*
         * System.out.println((int)Math.log10(121)+1);
         * 
         * System.out.println((int)Math.floor(121.0) % 10 );
         */

        /*
         * List<Integer> aux = new ArrayList<>();
         * 
         * aux.add(1);
         * aux.add(2);
         * aux.add(1);
         * List<Integer> aux2 = new ArrayList<>();
         * aux2.add(1);
         * aux2.add(2);
         * aux2.add(1);
         */
        // System.out.println(aux.hashCode() == aux2.hashCode());
        //
        // System.out.println(isPalindrome(12221));
        /*
         * System.out.println(Integer.MAX_VALUE+" "+Integer.MIN_VALUE);
         * 
         * System.out.println(myAtoi("-42"));
         * System.out.println(myAtoi("-0000000000042"));
         * System.out.println(myAtoi("                                   -042"));
         * System.out.println(myAtoi("42"));
         * System.out.println(myAtoi("133c111"));
         * System.out.println(myAtoi("                         0-1"));
         * System.out.println(myAtoi("-+12"));
         * System.out.println(myAtoi("+-12"));
         * System.out.println(myAtoi("+12"));
         * System.out.println(myAtoi("-"));
         * System.out.println(myAtoi("+1"));
         * System.out.println(myAtoi("    +0a32"));
         * System.out.println(myAtoi("words and 987"));
         * System.out.println(myAtoi("20000000000000000000"));
         * System.out.println(myAtoi(".1"));
         * System.out.println(myAtoi("  0000000000012345678"));
         */

        /*
         * printNumbers();
         * System.out.println(zigZagConversion("ABCD", 2));
         */
        /*
         * int[] nums = {1,2,2,3};
         * 
         * System.out.println(subsetsWithDup(nums));
         */
        /*
         * int[] nums = {6,3,2,7,4,-1};
         * 
         * 
         * System.out.println(permute(nums));
         */

        /*
         * int[] arrs = {16,51,36,29,84,80,46,97,84,16};
         * int target = 171;
         * System.out.println(threeSumMulti(arrs, target));
         */

        /*
         * int[] targets = {2,3,1,2,4,3};
         * int target = 7;
         */
        /*
         * int[] targets = {1,1,1,1,1,1,1,1};
         * int target = 11;
         * 
         * System.out.println(minSubArrayLen(target, targets));
         */
        /*
         * int[] numss ={-7,-8,7,5,7,1,6,0};
         * int k = 4;
         */

        /*
         * int[] numss =
         * {7157,9172,7262,-9146,3087,5117,4046,7726,-1071,6011,5444,-48,-1385,-7328,
         * 3255,1600,586,-5160,-371,-5978,9837,3255,-6137,8587,-3403,9775,260,6016,9797,
         * 3371,2395,6851,2349,-7019,9318,1211,-3110,8735,-7507,1784,7400,-5799,3169,-
         * 7696,-8991,-2222,-9434,-4490,4034,-831,-9656,5488,-4395,9339,4104,-9058,-4072
         * ,-1172,1758,6878,-5570,-6380,9550,-9389,1411,2298,3516,551,9196,5215,-237,-
         * 4146,1682,4418,-4639,7759,9593,-9588,3041,9208,-7331,-797,-2529,7738,-2944,
         * 4351,5091,-9448,-5404,6200,-1425,-3983,678,8456,-8085,5162,7165,4692,-494,-
         * 9249,8514,521,-8835,6745,-5775,-575,1876,-5464,5053,5567,3456,5873,1965,4316,
         * 2126,9462,-59,6544,-1547,7015,-8928,-3903,-3020,5865,-9479,6723,9214,5705,
         * 5136,7725,945,-1995,-2288,4579,7103,9938,4495,-730,-3180,7717,6824,794,-894,-
         * 1439,-1641,-4577,9362,-8817,-6035,-7980,-1278,-1928,-5390,-2342,1189,-2340,
         * 4788,-1814,5927,3115,9017,6801,7884,-5719,5992,7477,-486,-2734,-1557,3169,
         * 5288,-8295,-5651,2491,-3394,8302,-8822,5638,7654,7350,9884,-5392,881,-4874,
         * 5582,8309,-8514,2682,-6081,5602,4963,3538,9558,-6401,-2641,6223,-7107,-2772,
         * 5873,78,-7934,-7641,7872,7901,7436,-3815,-1540,-3387,3558,-8030,-6637,9609,
         * 8594,83,7984,-3286,7211,5877,-8655,6700,9855,-7521,903,1024,4051,4044,4044,
         * 8650,-2932,-134,-8167,-5338,-1014,391,1913,-9914,-9100,7108,-9250,1705,5615,
         * 641,6809,6619,7782,9062,3030,603,-2528,-5493,-1237,8428,1231,6701,202,641,-
         * 5351,-5366,-3347,7659,-3953,5518,1575,-3514,999,-6631,-934,-1119,1749,-9533,-
         * 8528,-9425,-9138,-6498,-1546,-8501,7668,-8135,-6234,7236,1722,-7690,7339,-
         * 5205,698,3680,7741,-9067,8739,-7658,-2518,3967,-128,620,-4571,780,5989,-6220,
         * -1932,6629,-733,-6978,-68,-3295,9075,-297,7648,-7645,2301,-4641,-8443,6935,-
         * 6257,7067,-9046,5474,6833,6924,8516,-213,-9210,-9605,-5798,4710,-9258,-7736,
         * 944,5194,-7465,5978,-6840,3735,4392,9218,-5571,2944,-5864,2995,-5234,5036,-
         * 4999,-9883,5493,4646,9574,3528,291,-4799,-3099,7639,5144,-2560,-7573,433,2464
         * ,-3484,4673,3283,-6459,-1194,8122,7314,-3389,-1899,8362,-1046,-1751,-2140,
         * 7642,-6274,-8056,3925,-397,1641,5762,8099,-9683,2533,1333,3295,7413,-8538,-
         * 8585,8412,1958,-8487,7248,4987,-6079,9427,-6207,-7873,688,224,6792,-4150,3345
         * ,826,1885,6463,-5269,3068,9649,-1354,3159,4975,514,-3071,-4355,-1615,9427,
         * 8343,978,7914,-1876,1160,-898,-8431,6245,8760,8514,9857,9505,-3602,-4124,-
         * 4124,209,855,-253,-7232,-7598,6813,-565,-8739,2886,3289,-4339,7846,-3820,3001
         * ,-3235,-3146,-2535,-1444,8976,-8434,8190,-4185,5847,-1020,-6020,-3935,-4267,
         * 2030,6882,-7707,-5213,5284,-2061,-325,2911,2346,1080,-2111,-4929,-9101,1548,-
         * 4817,-7526,2688,-3589,-4414,6269,-1423,-6735,-7204,-6624,-7561,7775,-2650,-
         * 6843,735,3824,4592,-5199,-1922,1757,5662,-1272,4208,400,2883,720,9179,1056,
         * 3310,-7095,-3834,-2683,4422,-2599,-6124,1449,-5001,-5874,-7396,9158,2926,4281
         * ,-9423,8492,-1542,1197,6023,-9627,4970,28,7002,5204,5292,3901,4640,2994,-4487
         * ,-2102,-4481,-5347,1164,6773,6277,5759,-4250,-3920,4843,7763,-791,8478,-7750,
         * 7243,-4640,6252,8699,2001,9799,-5555,-3183,-6124,4787,1378,-4618,3349,-5561,-
         * 2392,-1764,9774,-5698,1775,-9616,-6353,-3622,-4907,1356,5728,-1902,-3203,5268
         * ,4414,1096,-1268,-940,179,-7824,9845,6093,9096,-163,3713,-297,6100,6544,6167,
         * 6209,-5476,4519,6391,289,1823,7256,5528,9069,-4861,2571,-5339,2657,-1383,-
         * 3771,-4709,-1915,-8712,-816,2266,-8078,-2451,-6189,-5910,-8027,4915,-5900,-
         * 2979,2028,4015,-2885,8665,3121,8692,-2479,-2824,-5047,-3116,-5621,-7248,-1462
         * ,1114,-907,5481,6605,8767,-506,3412,-7848,7333,-634,3219,-3273,3031,-1867,
         * 1765,1522,-7747,-7195,-9110,6320,-3756,5207,1190,6370,-3143,6745,-2833,1926,-
         * 985,-3126,-9019,9744,-9202,8817,-3722,-2002,8111,4457,4973,4275,7125,3828,-
         * 3613,-3104,6544,6764,6585,-4240,-3961,-2756,-5445,-1143,-9788,-6964,3690,-
         * 1158,-6795,9726,7048,8414,-4774,8405,-8837,3163,-9265,877,-6371,-5901,5427,
         * 243,-8247,-2653,-2356,-1228,-3403,-9628,4430,1937,-8435,3876,-9615,-1366,-
         * 8793,2136,496,3957,-1316,822,7134,-8320,-8789,-33,1803,-2617,4625,-4334,-46,
         * 6870,-9895,-3381,-6536,7742,6356,-1725,2283,-2267,532,-3571,4288,-40,4714,
         * 2145,-8173,-9782,-2821,8418,7097,-7187,-2945,830,-1110,7886,-821,-3453,-4313,
         * -7945,7020,-2473,-4510,4867,-1992,3770,1031,6714,9721,-1399,-5297,-3545,-767,
         * -2432,-8088,-6801,1689,7271,673,9178,7565,8263,-213,6693,843,940,9793,7536,-
         * 1742,266,9280,-402,8335,5091,-3019,-3904,-6956,-7393,1053,9830,-403,6191,7652
         * ,-5990,-7726,741,-7996,-3664,-5601,9598,6603,3714,8336,5228,-3757,7069,-371,-
         * 9984,2625,-5485,-14,8394,7757,4705,-5743,-3141,6589,8246,7689,5709,9201,9740,
         * -5969,-3092,-5806,-1012,-7508,-9508,-9229,-6246,-5063,-8889,-4678,-7761,-4711
         * ,3076,-2699,224};
         * int k = 45;
         */
        // pruebaski();
        /*
         * int[] numss = { 1, 3, -1, -3, 5, 3, 6, 7 };
         * int k = 3;
         */
        /*
         * int[] numss={5,3,4};
         * int k = 1;
         */

        /*
         * int[]
         * numss={-1202,7068,-2705,3159,-398,786,5303,-9662,-109,-5256,4650,8345,2669,
         * 1465,-3552,5347,181,4830,-1018,-8237,-4305,4968,-1000,6762,-6620,-3714,6598,
         * 5681,4205,-4229,3879,-7038,-843,-9542,-8565,7175,-9772,-6923,-7681,46,-465,
         * 7551,-4129,7320,3035,9862,4845,-2688,3629,39,-3393,-1293,-3529,4773,-1466,-
         * 5074,2726,-8191,3948,3314,9705,3260,-1051,-3748,-2317,-5500,-2683,3465,3661,
         * 282,811,9982,-719,-5266,8031,6367,-5476,8762,-8380,-4587,7319,-8213,5100,5772
         * ,1470,7568,-4502,5221,9943,7754,-5680,-1315,516,-4648,-5349,-7531,-6459,-1804
         * ,9236,-1488,8296,216,6042,3655,9727,-7421,-9044,1754,-5846,3562,-8774,-2010,
         * 6587,707,-4524,-6424,-9977,4543,4414,-5942,-7926,7307,7430,4989,-836,-7486,-
         * 1114,2855,-1585,-4529,6097,-5446,-8334,4760,3856,5163,-3403,-4660,8997,5641,
         * 9133,-1149,8819,4569,8693,8762,-3662,-8779,-7500,-3206,6455,2988,-1719,3222,
         * 5285,-5353,2218,9778,2819,5020,-633,7670,1266,6912,-2276,-6220,-2925,-9810,
         * 1467,9826,-8802,9568,-3481,-6534,837,5802,1070,-857,2531,-3954,-754,-7219,
         * 7482,9506,8888,9469,-1574,9277,4272,-5004,1097,-9125,4043,341,-9764,-2048,
         * 7073,9455,2895,4680,2370,-9070,-5341,-3791,-9345,-9618,-1303,-6281,-9278,8299
         * ,1348,-3324,-6287,8402,675,-4075,-3545,7569,-7513,-4064,4181,4928,-4359,644,-
         * 4737,-562,9484,-4236,4165,-7783,8386,-1507,-5665,-2356,-9711,7275,6168,-7117,
         * -8124,-4373,538,8262,-2917,-3079,-5556,-7162,-9130,8760,-2914,6156,-886,-5288
         * ,-2986,870,-4815,2489,-6597,-864,-7425,-7773,-1538,2205,3488,8768,-7616,6045,
         * -5551,-4589,-4626,724,-4207,4682,1240,-3727,-7406,8433,-8006,-8774,3143,-8432
         * ,1955,-1840,3123,5899,6753,-9310,7266,-2150,-4088,3142,6370,-6001,-4766,-8007
         * ,4907,4385,9601,-3925,5503,394,7309,4318,-1713,-8997,-2839,1848,-787,-2743,-
         * 754,3443,8643,35,-3224,-5472,435,3916,-7864,-7392,-8496,3766,-3135,3300,-1613
         * ,-4373,-4375,-2045,763,8438,6930,-4186,-422,5766,-8255,-2209,-7796,-7047,8268
         * ,779,134,-7893,2870,-1691,-9568,-6078,-6125,4140,-757,-6069,-6905,1387,1007,
         * 9815,6128,5344,-1897,-1513,8967,145,-8649,-4234,7846,-1068,-4221,8569,3670,
         * 399,7222,3494,-4985,-484,7702,7159,-2844,6639,9702,-3468,-4091,5692,9112,7096
         * ,-6568,-9162,-2402,-8473,8428,5392,1244,-994,681,-3633,-2748,5203,2522,1290,-
         * 5566,8015,1593,-7248,-173,3344,7511,-604,4612,-9024,-1178,6610,7438,5786,7186
         * ,9147,-3525,3620,-2263,-2181,-5966,3410,5850,2103,-7645,2149,-4330,2668,4596,
         * -1189,2359,1587,-9086,-4473,9186,270,-7960,4856,4617,4374,51,7643,7946,2473,
         * 3813,-1783,6502,3871,-1440,1434,1650,1000,3266,-5874,-9296,-3093,4346,-6076,
         * 3663,959,5623,-8362,-9988,-6956,3257,-9702,-4119,8860,2932,-1053,-1354,1898,
         * 4736,4221,-8782,-1138,8141,2975,-8230,-4870,-1510,5618,-6556,-3505,-1477,9571
         * ,9814,-8088,7839,8371,-9265,9392,8364,-1385,-3676,-667,7562,2564,724,-5516,-
         * 9959,-8102,-1970,-3000,2295,9405,5965,-3218,6010,1067,2734,5528,-366,7328,
         * 4594,7320,-7491,6585,8741,1166,-8139,-6656,5510,8384,-2466,-9656,3767,2437,
         * 8638,2524,-8214,-1695,-5516,-9900,-9355,9930,2875,-9856,-5042,-4444,3574,5524
         * ,2583,1107,-2801,6824,-6156,5083,663,2178,-4962,-1961,5487,1199,2539,-9760,-
         * 2367,-5553,9240,-1283,9836,-3026,-1409,-3270,-2944,-3356,8075,968,7327,-7608,
         * -9929,3026,-8432,-8542,-7514,9644,-7411,5508,7608,-212,-1221,-9371,7944,-8788
         * ,-190,3278,7534,6508,5151,-2233,-5612,8522,5307,9136,4573,1468,1264,-5834,
         * 7043,1953,-5452,-6489,-4653,-5646,1787,-7525,-4227,-1763,7378,-6274,-1,9236,-
         * 9607,7803,4881,-7236,-4863,2634,-6935,-5012,-7445,2037,9333,-4374,6435,-2627,
         * 7397,-7166,9533,7799,-9842,388,-5030,31,-4047,-8042,626,2911,-7579,1945,853,-
         * 9632,-8468,-343,-957,8046,-5467,-6910,-9133,587,4738,2502,5115,-151,7323,-664
         * ,3205,8770,9584,-1615,5427,5714,7368,-2855,624,3432,-9769,-8064,-649,2346,398
         * ,9787,5232,-658,-2183,7548,-3126,5613,-6398,-3013,-1063,5968,-8088,2376,2275,
         * 7770};
         * int k = 456;
         */
        /*
         * int[] numss={1,3,1,2,0,5};
         * int k = 3;
         */

        /*
         * int[] numss={7,2,4};
         * int k = 2;
         */
        /*
         * int[] numss = {1,100,1,1,1,1,1,1,1,1,500};
         * int k=3;
         */

        /*
         * int[] numss =
         * {7157,9172,7262,-9146,3087,5117,4046,7726,-1071,6011,5444,-48,-1385,-7328,
         * 3255,1600,586,-5160,-371,-5978,9837,3255,-6137,8587,-3403,9775,260,6016,9797,
         * 3371,2395,6851,2349,-7019,9318,1211,-3110,8735,-7507,1784,7400,-5799,3169,-
         * 7696,-8991,-2222,-9434,-4490,4034,-831,-9656,5488,-4395,9339,4104,-9058,-4072
         * ,-1172,1758,6878,-5570,-6380,9550,-9389,1411,2298,3516,551,9196,5215,-237,-
         * 4146,1682,4418,-4639,7759,9593,-9588,3041,9208,-7331,-797,-2529,7738,-2944,
         * 4351,5091,-9448,-5404,6200,-1425,-3983,678,8456,-8085,5162,7165,4692,-494,-
         * 9249,8514,521,-8835,6745,-5775,-575,1876,-5464,5053,5567,3456,5873,1965,4316,
         * 2126,9462,-59,6544,-1547,7015,-8928,-3903,-3020,5865,-9479,6723,9214,5705,
         * 5136,7725,945,-1995,-2288,4579,7103,9938,4495,-730,-3180,7717,6824,794,-894,-
         * 1439,-1641,-4577,9362,-8817,-6035,-7980,-1278,-1928,-5390,-2342,1189,-2340,
         * 4788,-1814,5927,3115,9017,6801,7884,-5719,5992,7477,-486,-2734,-1557,3169,
         * 5288,-8295,-5651,2491,-3394,8302,-8822,5638,7654,7350,9884,-5392,881,-4874,
         * 5582,8309,-8514,2682,-6081,5602,4963,3538,9558,-6401,-2641,6223,-7107,-2772,
         * 5873,78,-7934,-7641,7872,7901,7436,-3815,-1540,-3387,3558,-8030,-6637,9609,
         * 8594,83,7984,-3286,7211,5877,-8655,6700,9855,-7521,903,1024,4051,4044,4044,
         * 8650,-2932,-134,-8167,-5338,-1014,391,1913,-9914,-9100,7108,-9250,1705,5615,
         * 641,6809,6619,7782,9062,3030,603,-2528,-5493,-1237,8428,1231,6701,202,641,-
         * 5351,-5366,-3347,7659,-3953,5518,1575,-3514,999,-6631,-934,-1119,1749,-9533,-
         * 8528,-9425,-9138,-6498,-1546,-8501,7668,-8135,-6234,7236,1722,-7690,7339,-
         * 5205,698,3680,7741,-9067,8739,-7658,-2518,3967,-128,620,-4571,780,5989,-6220,
         * -1932,6629,-733,-6978,-68,-3295,9075,-297,7648,-7645,2301,-4641,-8443,6935,-
         * 6257,7067,-9046,5474,6833,6924,8516,-213,-9210,-9605,-5798,4710,-9258,-7736,
         * 944,5194,-7465,5978,-6840,3735,4392,9218,-5571,2944,-5864,2995,-5234,5036,-
         * 4999,-9883,5493,4646,9574,3528,291,-4799,-3099,7639,5144,-2560,-7573,433,2464
         * ,-3484,4673,3283,-6459,-1194,8122,7314,-3389,-1899,8362,-1046,-1751,-2140,
         * 7642,-6274,-8056,3925,-397,1641,5762,8099,-9683,2533,1333,3295,7413,-8538,-
         * 8585,8412,1958,-8487,7248,4987,-6079,9427,-6207,-7873,688,224,6792,-4150,3345
         * ,826,1885,6463,-5269,3068,9649,-1354,3159,4975,514,-3071,-4355,-1615,9427,
         * 8343,978,7914,-1876,1160,-898,-8431,6245,8760,8514,9857,9505,-3602,-4124,-
         * 4124,209,855,-253,-7232,-7598,6813,-565,-8739,2886,3289,-4339,7846,-3820,3001
         * ,-3235,-3146,-2535,-1444,8976,-8434,8190,-4185,5847,-1020,-6020,-3935,-4267,
         * 2030,6882,-7707,-5213,5284,-2061,-325,2911,2346,1080,-2111,-4929,-9101,1548,-
         * 4817,-7526,2688,-3589,-4414,6269,-1423,-6735,-7204,-6624,-7561,7775,-2650,-
         * 6843,735,3824,4592,-5199,-1922,1757,5662,-1272,4208,400,2883,720,9179,1056,
         * 3310,-7095,-3834,-2683,4422,-2599,-6124,1449,-5001,-5874,-7396,9158,2926,4281
         * ,-9423,8492,-1542,1197,6023,-9627,4970,28,7002,5204,5292,3901,4640,2994,-4487
         * ,-2102,-4481,-5347,1164,6773,6277,5759,-4250,-3920,4843,7763,-791,8478,-7750,
         * 7243,-4640,6252,8699,2001,9799,-5555,-3183,-6124,4787,1378,-4618,3349,-5561,-
         * 2392,-1764,9774,-5698,1775,-9616,-6353,-3622,-4907,1356,5728,-1902,-3203,5268
         * ,4414,1096,-1268,-940,179,-7824,9845,6093,9096,-163,3713,-297,6100,6544,6167,
         * 6209,-5476,4519,6391,289,1823,7256,5528,9069,-4861,2571,-5339,2657,-1383,-
         * 3771,-4709,-1915,-8712,-816,2266,-8078,-2451,-6189,-5910,-8027,4915,-5900,-
         * 2979,2028,4015,-2885,8665,3121,8692,-2479,-2824,-5047,-3116,-5621,-7248,-1462
         * ,1114,-907,5481,6605,8767,-506,3412,-7848,7333,-634,3219,-3273,3031,-1867,
         * 1765,1522,-7747,-7195,-9110,6320,-3756,5207,1190,6370,-3143,6745,-2833,1926,-
         * 985,-3126,-9019,9744,-9202,8817,-3722,-2002,8111,4457,4973,4275,7125,3828,-
         * 3613,-3104,6544,6764,6585,-4240,-3961,-2756,-5445,-1143,-9788,-6964,3690,-
         * 1158,-6795,9726,7048,8414,-4774,8405,-8837,3163,-9265,877,-6371,-5901,5427,
         * 243,-8247,-2653,-2356,-1228,-3403,-9628,4430,1937,-8435,3876,-9615,-1366,-
         * 8793,2136,496,3957,-1316,822,7134,-8320,-8789,-33,1803,-2617,4625,-4334,-46,
         * 6870,-9895,-3381,-6536,7742,6356,-1725,2283,-2267,532,-3571,4288,-40,4714,
         * 2145,-8173,-9782,-2821,8418,7097,-7187,-2945,830,-1110,7886,-821,-3453,-4313,
         * -7945,7020,-2473,-4510,4867,-1992,3770,1031,6714,9721,-1399,-5297,-3545,-767,
         * -2432,-8088,-6801,1689,7271,673,9178,7565,8263,-213,6693,843,940,9793,7536,-
         * 1742,266,9280,-402,8335,5091,-3019,-3904,-6956,-7393,1053,9830,-403,6191,7652
         * ,-5990,-7726,741,-7996,-3664,-5601,9598,6603,3714,8336,5228,-3757,7069,-371,-
         * 9984,2625,-5485,-14,8394,7757,4705,-5743,-3141,6589,8246,7689,5709,9201,9740,
         * -5969,-3092,-5806,-1012,-7508,-9508,-9229,-6246,-5063,-8889,-4678,-7761,-4711
         * ,3076,-2699,224};
         * int k = 45;
         */

        /*
         * int[] numss = {6,2,3,500,400,4,2,6000,200,7,9999,9999,9999,9999,9999,9998};
         * int k = 3;
         */
        // System.out.println(maxSlidingWindow(numss, k));
        // target = 7, nums = [2,3,1,2,4,3]
        /*
         * int [] numsx = {2,3,1,2,4,3};
         * 
         * int targeting = 7;
         * 
         * System.out.println(minSubArrayLen(targeting, numsx));
         */
        // prueba();

        // char[][] cycles =
        // {{'a','a','a','a'},{'a','b','b','a'},{'a','b','b','a'},{'a','a','a','a'}};
        /*
         * char[][] cycles = { { 'c', 'c', 'c', 'a' }, { 'c', 'd', 'c', 'c' }, { 'c',
         * 'c', 'e', 'c' },
         * { 'f', 'c', 'c', 'c' } };
         */

        /*
         * char [][] cycles = {{'a','a','a'},
         * {'a','a','a'},
         * {'a','b','a'},
         * {'b','a','b'},
         * {'a','b','a'}};
         */

        /*
         * char[][] cycles = {{'c','a','d'},
         * {'a','a','a'},
         * {'a','a','d'},
         * {'a','c','d'},
         * {'a','b','c'}};
         */

        // System.out.println(containsCycle(cycles));
        // int [][] graphs={{0,10},{3,18},{5,5},{6,11},{11,14},{13,1},{15,1},{17,4}};
        // int[][] graphs = {{1,4},{2,4},{3,1},{3,2}}; //mal debe ser true

        // int[][] graphs = {{1,0},{2,1},{3,2},{1,3}};

        // int [][] graphs={{1,0},{0,2},{2,1}};
        // int[][] graphs =
        // {{1,0},{2,0},{2,1},{3,1},{3,2},{4,2},{4,3},{5,3},{5,4},{6,4},{6,5},{7,5},{7,6},{8,6},{8,7},{9,7},{9,8},{10,8},{10,9},{11,9},{11,10},{12,10},{12,11},{13,11},{13,12},{14,12},{14,13},{15,13},{15,14},{16,14},{16,15},{17,15},{17,16},{18,16},{18,17},{19,17},{19,18},{20,18},{20,19},{21,19},{21,20},{22,20},{22,21},{23,21},{23,22},{24,22},{24,23},{25,23},{25,24},{26,24},{26,25},{27,25},{27,26},{28,26},{28,27},{29,27},{29,28},{30,28},{30,29},{31,29},{31,30},{32,30},{32,31},{33,31},{33,32},{34,32},{34,33},{35,33},{35,34},{36,34},{36,35},{37,35},{37,36},{38,36},{38,37},{39,37},{39,38},{40,38},{40,39},{41,39},{41,40},{42,40},{42,41},{43,41},{43,42},{44,42},{44,43},{45,43},{45,44},{46,44},{46,45},{47,45},{47,46},{48,46},{48,47},{49,47},{49,48},{50,48},{50,49},{51,49},{51,50},{52,50},{52,51},{53,51},{53,52},{54,52},{54,53},{55,53},{55,54},{56,54},{56,55},{57,55},{57,56},{58,56},{58,57},{59,57},{59,58},{60,58},{60,59},{61,59},{61,60},{62,60},{62,61},{63,61},{63,62},{64,62},{64,63},{65,63},{65,64},{66,64},{66,65},{67,65},{67,66},{68,66},{68,67},{69,67},{69,68},{70,68},{70,69},{71,69},{71,70},{72,70},{72,71},{73,71},{73,72},{74,72},{74,73},{75,73},{75,74},{76,74},{76,75},{77,75},{77,76},{78,76},{78,77},{79,77},{79,78},{80,78},{80,79},{81,79},{81,80},{82,80},{82,81},{83,81},{83,82},{84,82},{84,83},{85,83},{85,84},{86,84},{86,85},{87,85},{87,86},{88,86},{88,87},{89,87},{89,88},{90,88},{90,89},{91,89},{91,90},{92,90},{92,91},{93,91},{93,92},{94,92},{94,93},{95,93},{95,94},{96,94},{96,95},{97,95},{97,96},{98,96},{98,97},{99,97}};
        // int[][] graphs = {{0,1},{2,3},{1,2},{3,0}};
        // int [][] graphs={{1,0},{0,1}};
        // System.out.println(canFinish(4, graphs));

        /*
         * System.out.println("Hello World!");
         * System.out.println("ESCALERA reves "+fib(3));
         */

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

        // String prueba = "ABC";
        // int numRows = 2;

        // System.out.println(convert(prueba, numRows));

        // int[] time = { 20, 40 };
        // System.out.println(numPairsDivisibleBy60(time));

        // System.out.println((60-(150%60)) %60);
        // System.out.println(150%60);

        /**
         * RETOMA ESTE PROBLERMA AL PENULTIMO O AL RESOLVER OTRO PROBLEMA
         * String[] words =
         * {"dd","aa","bb","dd","aa","dd","bb","dd","aa","cc","bb","cc","dd","cc"};
         * System.out.println(longestPalindrome(words));
         */

        // int[] caseD = { 10, 20, 40, 80 };
        // System.out.println(Math.absExact(-3 + 2) + "+");

        // System.out.println(canReorderDoubledV2(caseD));

        /**
         * int[] arr = { 7 , 1 , 5 , 3 , 6 , 4};
         * 
         * System.out.println(returnMax(arr));
         */

        // int[] testing = { 1, 2, 3, 1 };
        // System.out.println(Math.abs(0-3)); // es menor o igual a indexDiff que es
        // tres
        /*
         * char[][] arr = {{'.','.','4','.','.','.','6','3','.'},
         * {'.','.','.','.','.','.','.','.','.'},
         * {'5','.','.','.','.','.','.','9','.'},
         * {'.','.','.','5','6','.','.','.','.'},
         * {'4','.','3','.','.','.','.','.','1'},
         * {'.','.','.','7','.','.','.','.','.'},
         * {'.','.','.','5','.','.','.','.','.'},
         * {'.','.','.','.','.','.','.','.','.'},
         * {'.','.','.','.','.','.','.','.','.'}};
         */
        // System.out.println(containsNearbyAlmostDuplicate(testing,3,0));

        // System.out.println(arr[0].length + "<>" + arr.length);

        // System.out.println(isValidSudoku(arr));

        // String[] strs = {"eat","tea","tan","ate","nat","bat"};

        // System.out.println(groupAnagrams(strs));

        /*
         * String longest = "au";
         * System.out.println(lengthOfLongestSubstring(longest));
         */

        // String s2 ="abab",p2="ab";
        // System.out.println(findAnagrams(s, p));

        // int[] test = {1,1,2,2,3,3,4,4,5,5};

        // System.out.println(threeSumMulti(test, 8));

        /*
         * long start = System.nanoTime();
         * int[][] test = {{1,1,1},{1,2,2},{1,2,2}};
         * int [][] outputing = floodFill(test, 1, 1, 3);
         * long end = System.nanoTime();
         * 
         * for (int[] is : outputing) {
         * for (int is2 : is) {
         * System.out.print(is2);
         * }
         * System.out.println();
         * }
         * 
         * System.out.println((double) end-start/1000000);
         */

        /*
         * char[][] islands = {
         * {'1','1','1','1','0'},
         * {'1','1','0','1','0'},
         * {'0','0','0','0','0'},
         * {'0','0','0','0','0'}};
         */

        // System.out.println(numIslands(islands));

        /*
         * int[][] islandsC = {{0,0,1,0,0,0,0,1,0,0,0,0,0},
         * {0,0,0,0,0,0,0,1,1,1,0,0,0},
         * {0,1,1,0,1,0,0,0,0,0,0,0,0},
         * {0,1,0,0,1,1,0,0,1,0,1,0,0},
         * {0,1,0,0,1,1,0,0,1,1,1,0,0},
         * {0,0,0,0,0,0,0,0,0,0,1,0,0},
         * {0,0,0,0,0,0,0,1,1,1,0,0,0},
         * {0,0,0,0,0,0,0,1,1,0,0,0,0}};
         */

        // System.out.println(maxAreaOfIsland(islandsC));

        // int[] nums = { 1, 2, 2 };
        // System.out.println(subsetsWithDup(nums));

        // System.out.println(permute(new int[]{1,2}));

        // System.out.println("ESCALERA "+climbStairs(6));

    }
}

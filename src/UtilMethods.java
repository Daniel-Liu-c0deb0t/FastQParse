import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

public class UtilMethods {
	private static HashMap<Character, Character> complements = new HashMap<Character, Character>();
	
	static{
		//complements for each bp
		complements.put('A', 'T');
		complements.put('T', 'A');
		complements.put('C', 'G');
		complements.put('G', 'C');
		complements.put('N', 'N');
	}
	
	public static void defaultExceptionHandler(PrintWriter logWriter, Throwable e){
		if(logWriter != null){
			e.printStackTrace(logWriter);
			logWriter.flush();
		}else{
			e.printStackTrace();
		}
		System.exit(1);
	}
	
	//create a byte array from a quality score string
	public static byte[] qScoreToByteArray(String q){
		byte[] bytes = new byte[q.length()];
		for(int i = 0; i < q.length(); i++){
			bytes[i] = (byte)toQScore(q.charAt(i));
		}
		return bytes;
	}
	
	//convert a byte array to a quality score string
	public static String byteArrayToQScore(byte[] arr){
		char[] chars = new char[arr.length];
		for(int i = 0; i < arr.length; i++){
			chars[i] = toQualityChar(arr[i]);
		}
		return new String(chars);
	}
	
	//count number of N bp
	public static int countN(String s){
		int count = 0;
		for(int i = 0; i < s.length(); i++){
			if(Character.toUpperCase(s.charAt(i)) == 'N'){
				count++;
			}
		}
		return count;
	}
	
	//percentage of N bp
	public static double percentN(String s){
		if(s.length() == 0)
			return 0;
		return (double)countN(s) / (double)s.length();
	}
	
	//expected number of bases that are wrong
	//sum of all error percentages
	//skips characters from start if offset > 0
	public static double toError(String quality, int offset){
		double result = 0.0;
		for(int i = offset; i < quality.length(); i++){
			result += toError(quality.charAt(i));
		}
		return result;
	}
	
	//get error percentage from quality character
	public static double toError(char quality){
		return qScoreToError(toQScore(quality));
	}
	
	//average of quality scores
	//skip characters from start if offset > 0
	public static double toQScore(String quality, int offset){
		double result = 0.0;
		for(int i = offset; i < quality.length(); i++){
			result += toQScore(quality.charAt(i));
		}
		return result / quality.length();
	}
	
	//positive integer of quality, greater = better (less likely it is wrong)
	//quality of one character
	public static int toQScore(char quality){
		return quality - '!'; //! is 33rd ASCII
	}
	
	//converts quality score (int) to a error percentage (double)
	public static double qScoreToError(int qScore){
		return Math.pow(10.0, qScore / -10.0);
	}
	
	//get error percentage if the bases are the same
	public static double sameError(double error1, double error2){
		return (error1 * error2 / 3.0) / (1.0 - error1 - error2 + 4.0 * error1 * error2 / 3.0);
	}
	
	//get error percentage if the bases are different
	public static double differentError(double error1, double error2){
		double temp = error1;
		error1 = Math.min(error1, error2);
		error2 = Math.max(temp, error2);
		return error1 * (1 - error2 / 3.0) / (error1 + error2 - 4.0 * error1 * error2 / 3.0);
	}
	
	//get quality char from error percentage
	public static char toQualityChar(double error){
		return toQualityChar(errorToQScore(error));
	}
	
	//convert error percentage (double) to quality score (int)
	public static int errorToQScore(double error){
		return (int)Math.round(-10.0 * Math.log10(error));
	}
	
	//converts an integer into a quality character
	//0 converts to '!'
	public static char toQualityChar(int quality){
		if(quality + '!' >= '!' && quality + '!' <= '~'){ //check if the quality is in the bounds
			return (char)(quality + '!');
		}
		return (quality + '!' < '!') ? '!' : '~';
	}
	
	public static BitSet toBit(String s){
		BitSet set = new BitSet(s.length() * 3);
		for(int i = 0; i < s.length(); i++){
			if(Character.toUpperCase(s.charAt(i)) == 'A'){ //100
				set.set(i * 3);
			}else if(Character.toUpperCase(s.charAt(i)) == 'T'){ //010
				set.set(i * 3 + 1);
			}else if(Character.toUpperCase(s.charAt(i)) == 'C'){ //110
				set.set(i * 3);
				set.set(i * 3 + 1);
			}else if(Character.toUpperCase(s.charAt(i)) == 'G'){ //001
				set.set(i * 3 + 2);
			}else if(Character.toUpperCase(s.charAt(i)) == 'N'){ //101
				set.set(i * 3);
				set.set(i * 3 + 2);
			}
		}
		return set;
	}
	
	public static String toSequence(BitSet set){
		StringBuilder builder = new StringBuilder();
		for(int i = 0;; i++){
			if(set.get(i * 3) && !set.get(i * 3 + 1) && !set.get(i * 3 + 2)){ //100
				builder.append('A');
			}else if(!set.get(i * 3) && set.get(i * 3 + 1) && !set.get(i * 3 + 2)){ //010
				builder.append('T');
			}else if(set.get(i * 3) && set.get(i * 3 + 1) && !set.get(i * 3 + 2)){ //110
				builder.append('C');
			}else if(!set.get(i * 3) && !set.get(i * 3 + 1) && set.get(i * 3 + 2)){ //001
				builder.append('G');
			}else if(set.get(i * 3) && !set.get(i * 3 + 1) && set.get(i * 3 + 2)){ //101
				builder.append('N');
			}else{
				break;
			}
		}
		return builder.toString();
	}
	
	//recursive deletion of everything in a folder and the folder itself
	public static void deleteFolder(File folder){
		File[] files = folder.listFiles();
		if(files != null){
			for(File f : files){
				if(f.isDirectory()){
					deleteFolder(f);
				}else{
					f.delete();
				}
			}
		}
		folder.delete();
	}
	
	//make a string that consists of a bunch of one character
	public static String makeStr(char c, int n){
		char[] result = new char[n];
		Arrays.fill(result, c);
		return new String(result);
	}
	
	//fuzzy match where N can stand for any character
	public static boolean matchWithN(String a, String b, double max, boolean indel, boolean wildcard){
		return distWithN(a, b, indel, wildcard) <= (max < 0.0 ? (-max * b.length()) : max);
	}
	
	//distance between two strings
	//N can stand for any character
	//support for insertions, deletions, and substitutions
	public static int distWithN(String a, String b, boolean indel, boolean wildcard){
		if(a.isEmpty() || b.isEmpty())
			return Math.max(a.length(), b.length());
		
		if(indel){ //Wagner-Fischer algorithm
			if(a.length() > b.length()){ //use the smaller string length as the array size
				String temp = a;
				a = b;
				b = temp;
			}
			
			int[] curr = new int[a.length() + 1]; //two arrays to save space
			int[] prev = new int[a.length() + 1];
			
			for(int i = 1; i <= b.length(); i++){
				curr[0] = i;
				for(int j = 1; j <= a.length(); j++){
					if(Character.toUpperCase(b.charAt(i - 1)) == Character.toUpperCase(a.charAt(j - 1)) ||
							(wildcard && (Character.toUpperCase(b.charAt(i - 1)) == 'N' || Character.toUpperCase(a.charAt(j - 1)) == 'N'))){
						curr[j] = i == 1 ? (j - 1) : prev[j - 1]; //if bp equals or they equal N
					}else{ //insertion, deletion, and substitution
						curr[j] = Math.min(Math.min((i == 1 ? j : prev[j]) + 1, curr[j - 1] + 1), (i == 1 ? (j - 1) : prev[j - 1]) + 1);
					}
					if(i < b.length())
						prev[j - 1] = curr[j - 1]; //update prev array with previously calculated elements in curr array
				}
				if(i < b.length())
					prev[a.length()] = curr[a.length()];
			}
			
			return curr[a.length()];
		}else{ //count mismatches (substitution distance)
			int minLen = Math.min(a.length(), b.length());
			int wrong = Math.max(a.length(), b.length()) - minLen;
			
			for(int i = 0; i < minLen; i++){
				if(wildcard && (Character.toUpperCase(a.charAt(i)) == 'N' || Character.toUpperCase(b.charAt(i)) == 'N'))
					continue;
				if(Character.toUpperCase(a.charAt(i)) != Character.toUpperCase(b.charAt(i))){
					wrong++;
				}
			}
			
			return wrong;
		}
	}
	
	//fuzzy search for string b in a
	//returns a list of string end positions + 1
	//set offset to Integer.MAX_VALUE for searching
	//set minOverlap to Integer.MAX_VALUE to make sure that the match only appears within the string to be searched
	public static ArrayList<Match> searchWithN(String a, String b, double max, int offset, boolean indel, boolean bestOnly, int minOverlap, boolean wildcard){
		if(b.isEmpty())
			return new ArrayList<Match>(Arrays.asList(new Match(0, 0, 0)));
		
		//try hamming (substitution only) search first
		ArrayList<Match> hammingResults = searchWithNHamming(a, b, max, offset, bestOnly, minOverlap, wildcard);
		if(!hammingResults.isEmpty()){
			return hammingResults;
		}
		
		if(indel){ //Wagner-Fischer algorithm, with the ability to search and match different lengths
			int[][] curr = new int[b.length() + 1][3]; //{insertion, deletion, substitution}
			int[][] prev = new int[b.length() + 1][3];
			ArrayList<Match> result = new ArrayList<Match>();
			int min = Integer.MAX_VALUE;
			int end = Math.min((int)(max < 0.0 ? (-max * b.length()) : max) + b.length() - (minOverlap > b.length() ? b.length() : minOverlap) + 1, b.length());
			
//			for(int i = 0; i <= b.length(); i++)
//				System.out.print(Math.max(0, minOverlap < b.length() ? (i - b.length() + minOverlap) : i) + " ");
//			System.out.println();
			for(int i = 1; i <= a.length(); i++){
				curr[0] = new int[]{0, Math.max(0, i - 1 - offset), 0};
				for(int j = 1; j <= end; j++){
					if(Character.toUpperCase(b.charAt(j - 1)) == Character.toUpperCase(a.charAt(i - 1)) ||
							(wildcard && (Character.toUpperCase(b.charAt(j - 1)) == 'N' || Character.toUpperCase(a.charAt(i - 1)) == 'N'))){
						curr[j] = i == 1 ? new int[]{Math.max(0, minOverlap < b.length() ? (j - 1 - b.length() + minOverlap) : j - 1), 0, 0} : copy(prev[j - 1]);
					}else{
						int sub = i == 1 ? Math.max(0, minOverlap < b.length() ? (j - 1 - b.length() + minOverlap) : j - 1) : sum(prev[j - 1]);
						int ins = sum(curr[j - 1]);
						int del = i == 1 ? Math.max(0, minOverlap < b.length() ? (j - b.length() + minOverlap) : j) : sum(prev[j]);
						if(sub <= ins && sub <= del){
							curr[j] = i == 1 ? new int[]{Math.max(0, minOverlap < b.length() ? (j - 1 - b.length() + minOverlap) : j - 1), 0, 1} : new int[]{prev[j - 1][0], prev[j - 1][1], prev[j - 1][2] + 1};
						}else if(ins <= sub && ins <= del){
							curr[j] = new int[]{curr[j - 1][0] + 1, curr[j - 1][1], curr[j - 1][2]};
						}else if(del <= sub && del <= ins){
							curr[j] = i == 1 ? new int[]{Math.max(0, minOverlap < b.length() ? (j - b.length() + minOverlap) : j), 1, 0} : new int[]{prev[j][0], prev[j][1] + 1, prev[j][2]};
						}
					}
					prev[j - 1] = copy(curr[j - 1]);
//					System.out.print(sum(curr[j - 1]) + " ");
				}
				prev[end] = copy(curr[end]);
//				System.out.println(sum(curr[end]));
				while(end >= 0 && sum(curr[end]) > (max < 0.0 ? (-max * b.length()) : max)){
					end--;
				}
				if(end == b.length()){
					int index = i;
					int length;
					if(index < b.length()){
						length = index;
					}else{
						length = b.length();
					}
					length += curr[b.length()][0] - curr[b.length()][1];
					if(length >= minOverlap && sum(curr[b.length()]) <= (max < 0.0 ? (-max * length) : max)){
						if(!bestOnly || sum(curr[b.length()]) <= min){
							result.add(new Match(index, sum(curr[b.length()]), length));
							min = sum(curr[b.length()]);
						}
					}
				}else{
					end++;
				}
			}
			
			if(bestOnly && !result.isEmpty()){ //if only find best, then return only the matches with shortest distance
				ArrayList<Match> result2 = new ArrayList<Match>();
				for(int i = result.size() - 1; i >= 0; i--){
					if(result.get(i).edits == min){
						result2.add(result.get(i));
					}else{
						break;
					}
				}
				return result2;
			}
			return result;
		}else{
			return hammingResults;
		}
	}
	
	//very simple substitution only search
	public static ArrayList<Match> searchWithNHamming(String a, String b, double max, int offset, boolean bestOnly, int minOverlap, boolean wildcard){
		if(a.length() < b.length())
			return new ArrayList<Match>();
		
		ArrayList<Match> result = new ArrayList<Match>();
		int min = Integer.MAX_VALUE;
		
		for(int i = Math.min(minOverlap, b.length()); i <= a.length(); i++){
			int dist = distWithN(a.substring(Math.max(0, i - b.length()), i), b.substring(i < b.length() ? (b.length() - i) : 0), false, wildcard);
			int index = i;
			int length;
			if(index < b.length())
				length = index;
			else
				length = b.length();
			if(dist <= (max < 0.0 ? (-max * length) : max)){
				if(!bestOnly || dist <= min){
					result.add(new Match(index, dist, length));
					min = dist;
				}
			}
		}
		
		if(bestOnly && !result.isEmpty()){
			ArrayList<Match> result2 = new ArrayList<Match>();
			for(int i = result.size() - 1; i >= 0; i--){
				if(result.get(i).edits == min){
					result2.add(result.get(i));
				}else{
					break;
				}
			}
			return result2;
		}
		return result;
	}
	
	public static int sum(int[] arr){
		int sum = 0;
		for(int i = 0; i < arr.length; i++){
			sum += arr[i];
		}
		return sum;
	}
	
	public static int[] copy(int[] arr){
		int[] result = new int[arr.length];
		for(int i = 0; i < arr.length; i++){
			result[i] = arr[i];
		}
		return result;
	}
	
	//format elapsed time like 00:00:00.000 (hours, minutes, seconds, milliseconds)
	public static String formatElapsedTime(long time){
		return String.format("%02d:%02d:%02d.%03d", time / (3600 * 1000), time / (60 * 1000) % 60, time / 1000 % 60, time % 1000);
	}
	
	//current memory usage in GB as a double
	public static double currentMemoryUsage(){
		return (double)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1024.0 / 1024.0 / 1024.0;
	}
	
	//reverse a string
	public static String reverse(String s){
		char[] arr = s.toCharArray();
		char temp;
		
		for(int i = 0; i < arr.length / 2; i++){
			temp = arr[i];
			arr[i] = arr[arr.length - 1 - i];
			arr[arr.length - 1 - i] = temp;
		}
		return new String(arr);
	}
	
	//find the complement of a string of A, T, C, G, and N
	public static String complement(String s){
		char[] result = new char[s.length()];
		
		for(int i = 0; i < s.length(); i++){
			result[i] = complement(s.charAt(i));
		}
		return new String(result);
	}
	
	//find the complement of one character that is either A, T, C, G, or N
	public static char complement(char c){
		if(Character.isUpperCase(c))
			return complements.get(c);
		else
			return Character.toLowerCase(complements.get(Character.toUpperCase(c)));
	}
	
	//find the reverse complement of a string of A, T, C, G, and N
	public static String reverseComplement(String s){
		char[] arr = s.toCharArray();
		char temp;
		
		for(int i = 0; i < arr.length / 2; i++){
			temp = arr[i];
			arr[i] = complement(arr[arr.length - 1 - i]);
			arr[arr.length - 1 - i] = complement(temp);
		}
		if(arr.length % 2 != 0)
			arr[arr.length / 2] = complement(arr[arr.length / 2]);
		
		return new String(arr);
	}
	
	//merge two reads
	//increase quality if two base pairs are equal
	//decrease quality if two base pairs are not equal
	public static String[] mergeReads(String s1, String q1, String s2, String q2, double editMax, boolean algorithm, boolean wildcard){
		s2 = reverseComplement(s2);
		q2 = reverse(q2);
		
		int start = -1;
		
		if(algorithm){ //quality based matching
			double maxScore = 0;
			double e1;
			double e2;
			for(int i = 0; i < s1.length(); i++){
				double score = 0;
				for(int j = 0; j < Math.min(s2.length(), s1.length() - i); j++){
					e1 = toError(q1.charAt(i + j));
					e2 = toError(q2.charAt(j));
					
					if(wildcard && (Character.toUpperCase(s1.charAt(i + j)) == 'N' || Character.toUpperCase(s2.charAt(j)) == 'N')){
						score += 0.5;
					}else if(s1.charAt(i + j) == s2.charAt(j)){
						score += (1.0 - e1) * (1.0 - e2) + e1 * e2 / 3.0;
					}else{
						score -= 1.0 - ((1.0 - e2) * e1 / 3.0 + (1.0 - e1) * e2 / 3.0 + e1 * e2 / 2.0);
					}
				}
				if(score > maxScore){
					start = i;
					maxScore = score;
				}
			}
		}else{ //find match with least edit distance
			int minError = Integer.MAX_VALUE;
			int error = 0;
			int maxNonError = 0;
			int nonError = 0;
			for(int i = 0; i < s1.length(); i++){
				error = distWithN(s1.substring(i, Math.min(s1.length(), i + s2.length())), s2.substring(0, Math.min(s2.length(), s1.length() - i)), false, wildcard);
				nonError = Math.min(s2.length(), s1.length() - i) - error;
				if(error <= (editMax < 0.0 ? (-editMax * Math.min(s2.length(), s1.length() - i)) : editMax) && nonError >= maxNonError && (nonError > maxNonError || error < minError)){
					start = i;
					maxNonError = nonError;
					minError = error;
				}
			}
		}
		
		if(start == -1){
			return new String[]{s1 + s2, q1 + q2};
		}else{ //reconstruct the new, merged string
			StringBuilder sb1 = new StringBuilder(); //DNA data
			StringBuilder sb2 = new StringBuilder(); //quality
			for(int i = 0; i < Math.max(s1.length(), start + s2.length()); i++){
				if(i < start){
					sb1.append(s1.charAt(i));
					sb2.append(q1.charAt(i));
				}else if(i < Math.min(s1.length(), start + s2.length())){
					if((wildcard && (Character.toUpperCase(s1.charAt(i)) == 'N' || Character.toUpperCase(s2.charAt(i - start)) == 'N')) || s1.charAt(i) == s2.charAt(i - start)){
						sb1.append(s1.charAt(i));
						sb2.append(toQualityChar(sameError(toError(q1.charAt(i)), toError(q2.charAt(i - start)))));
					}else{
						sb1.append(toError(q1.charAt(i)) < toError(q2.charAt(i - start)) ? s1.charAt(i) : s2.charAt(i - start));
						sb2.append(toQualityChar(differentError(toError(q1.charAt(i)), toError(q2.charAt(i - start)))));
					}
				}else if(s1.length() < start + s2.length()){
					sb1.append(s2.charAt(i - start));
					sb2.append(q2.charAt(i - start));
				}else{
					sb1.append(s1.charAt(i));
					sb2.append(q1.charAt(i));
				}
			}
			return new String[]{sb1.toString(), sb2.toString()};
		}
	}
	
	//remove adapters from sequence and quality strings
	public static String[] removeAdapters(String s, String q, ArrayList<Adapter> adapters, double editMax, int minOverlap, int maxOffset, boolean indel, boolean mode, boolean wildcard){
		if(mode){ //allow the adapter to hang off each end of the read (deprecated)
			for(int i = 0; i < adapters.size(); i++){
				Adapter a = adapters.get(i);
				if(a.anchored){
					if(a.isStart){
						if(matchWithN(s.substring(0, a.str.length()), a.str, editMax, indel, wildcard)){
							s = s.substring(a.str.length());
							q = q.substring(a.str.length());
						}
					}else{
						if(matchWithN(s.substring(s.length() - a.str.length()), a.str, editMax, indel, wildcard)){
							s = s.substring(0, s.length() - a.str.length());
							q = q.substring(0, q.length() - a.str.length());
						}
					}
				}else{
					int maxNonError = 0;
					int minError = Integer.MAX_VALUE;
					int nonError;
					int error = 0;
					int bestIndex = -1;
					for(int j = a.isStart ? minOverlap : s.length() + a.str.length() - minOverlap; a.isStart ? j <= s.length() + a.str.length() - minOverlap : j >= minOverlap; j += a.isStart ? 1 : -1){
						if(a.isStart && j <= a.str.length()){
							error = distWithN(s.substring(0, j), a.str.substring(a.str.length() - j), indel, wildcard);
							nonError = j - error;
							if(error <= (editMax < 0.0 ? (-editMax * j) : editMax) && nonError >= maxNonError && (nonError > maxNonError || error < minError)){
								bestIndex = j;
								maxNonError = nonError;
								minError = error;
							}
						}else if(j <= s.length()){
							error = distWithN(s.substring(j - Math.min(j, a.str.length()), j), a.str, indel, wildcard);
							nonError = Math.min(j, a.str.length()) - error;
							if(error <= (editMax < 0.0 ? (-editMax * Math.min(j, a.str.length())) : editMax) && nonError >= maxNonError && (nonError > maxNonError || error < minError)){
								bestIndex = j;
								maxNonError = nonError;
								minError = error;
							}
						}else if(!a.isStart){
							error = distWithN(s.substring(j - a.str.length()), a.str.substring(0, a.str.length() - j + s.length()), indel, wildcard);
							nonError = a.str.length() - j + s.length() - error;
							if(error <= (editMax < 0.0 ? (-editMax * (a.str.length() - j + s.length())) : editMax) && nonError >= maxNonError && (nonError > maxNonError || error < minError)){
								bestIndex = j;
								maxNonError = nonError;
								minError = error;
							}
						}
					}
					
					if(bestIndex != -1){
						if(a.isStart){
							s = s.substring(bestIndex);
							q = q.substring(bestIndex);
						}else{
							s = s.substring(0, bestIndex - a.str.length());
							q = q.substring(0, bestIndex - a.str.length());
						}
					}
				}
			}
			
			return new String[]{s, q};
		}else{ //search with flexible starting and ending location
			for(int i = 0; i < adapters.size(); i++){
				Adapter a = adapters.get(i);
				
				ArrayList<Match> matches;
				if(a.anchored){
					matches = searchWithN(a.isStart ? s.substring(0, Math.min(a.str.length() + (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax), s.length())) :
						reverse(s).substring(0, Math.min(a.str.length() + (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax), s.length())), a.isStart ? a.str : reverse(a.str), editMax, 0, indel, true, Integer.MAX_VALUE, wildcard);
				}else{ //because searchWithN can only find ending locations, the 3' adapters need to be reversed along with the read
					ArrayList<Match> tempMatches = searchWithN(a.isStart ? s.substring(0, Math.min(maxOffset + a.str.length() + (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax), s.length())) :
						reverse(s).substring(0, Math.min(maxOffset + a.str.length() + (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax), s.length())), a.isStart ? a.str : reverse(a.str), editMax, maxOffset, indel, false, minOverlap, wildcard);
					matches = new ArrayList<Match>();
					
					int minEdit = Integer.MAX_VALUE;
					int maxLength = 0;
					int matchIndex = -1;
					for(int j = 0; j < tempMatches.size(); j++){
						Match m = tempMatches.get(j);
						if(m.correctLength >= maxLength && (m.correctLength > maxLength || m.edits < minEdit)){
							matchIndex = j;
							maxLength = m.correctLength;
							minEdit = m.edits;
						}
					}
					if(matchIndex != -1){
						matches.add(tempMatches.get(matchIndex));
					}
				}
				if(!matches.isEmpty()){
					if(a.isStart){
						s = s.substring(matches.get(matches.size() - 1).end);
						q = q.substring(matches.get(matches.size() - 1).end);
					}else{
						s = s.substring(0, s.length() - matches.get(matches.size() - 1).end); //reverse the index to get the correct index
						q = q.substring(0, q.length() - matches.get(matches.size() - 1).end);
					}
				}
			}
			
			return new String[]{s, q};
		}
	}
	
	//quality trim method 1
	public static String[] qualityTrim1(String s, String q, int minQuality, boolean trimLeft, int length){
		if(minQuality == 0){
			return new String[]{s, q};
		}
		
		for(int i = trimLeft ? 0 : s.length() - 1; trimLeft ? i < s.length() : i >= 0; i += trimLeft ? 1 : -1){
			int count = 0;
			double avg = 0.0;
			int q1 = toQScore(q.charAt(i));
			for(int j = (length == Integer.MAX_VALUE ? 0 : Math.max(i - length, 0)); j < (length == Integer.MAX_VALUE ? s.length() : Math.min(i + length + 1, s.length())); j++){ //average nearby quality scores
				int q2 = toQScore(q.charAt(j));
				if(q2 <= q1){ //only average with nearby bp with worse quality
					count++;
					avg += q2;
				}
			}
			if(avg / count >= minQuality){ //trimming is complete as soon as the average exceeds the threshold
				if(trimLeft){
					return new String[]{s.substring(i), q.substring(i)};
				}else{
					return new String[]{s.substring(0, i + 1), q.substring(0, i + 1)};
				}
			}
		}
		return new String[]{"", ""};
	}
	
	//quality trim method 2
	public static String[] qualityTrim2(String s, String q, int minQuality, boolean trimLeft, int length){
		if(minQuality == 0){
			return new String[]{s, q};
		}
		
		long count = 0l;
		for(int i = trimLeft ? 0 : s.length() - 1; trimLeft ? i < s.length() : i >= 0; i += trimLeft ? 1 : -1){
			//sliding window quality trim
			//the current count is previous count - the last element + the next element
			if(trimLeft) //subtract last element
				count -= (length == Integer.MAX_VALUE || (i - length) < 0) ? 0 : (toQScore(q.charAt(i - length)) - minQuality);
			else //subtract last element
				count -= (length == Integer.MAX_VALUE || (i + length) >= s.length()) ? 0 : (toQScore(q.charAt(i + length)) - minQuality);
			count += toQScore(q.charAt(i)) - minQuality; //add next element
			if(count >= 0l){
				if(trimLeft){
					return new String[]{s.substring(i), q.substring(i)};
				}else{
					return new String[]{s.substring(0, i + 1), q.substring(0, i + 1)};
				}
			}
		}
		return new String[]{"", ""};
	}
	
	//trim N based on a percentage from both sides of a read
	public static String[] trimN(String s, String q, double maxPercent){
		if(maxPercent > 1.0)
			return new String[]{s, q};
		
		int bestIndex = -1;
		int count = 0;
		double bestPercent = 0.0;
		boolean prevN = true;
		for(int i = 0; i < s.length(); i++){
			if(Character.toUpperCase(s.charAt(i)) == 'N'){
				if(prevN) //skip initial N
					continue;
				count++; //if not initial N, then count it
			}else{
				if(prevN){ //first bp that is not N
					bestIndex = i;
					prevN = false;
				}
			}
			double percent = (double)count / (i + 1.0); //percentage of N from the start of the read to current index
			if(percent >= maxPercent && percent >= bestPercent){
				bestIndex = i + 1;
				bestPercent = percent;
			}
		}
		
		if(prevN) //if read is all N
			return new String[]{"", ""};
		
		if(bestIndex != -1){
			s = s.substring(bestIndex);
			q = q.substring(bestIndex);
		}
		
		//reverse of the above
		bestIndex = -1;
		count = 0;
		bestPercent = 0.0;
		prevN = true;
		for(int i = s.length() - 1; i >= 0; i--){
			if(Character.toUpperCase(s.charAt(i)) == 'N'){
				if(prevN)
					continue;
				count++;
			}else{
				if(prevN){
					bestIndex = i;
					prevN = false;
				}
			}
			double percent = (double)count / (double)(s.length() - i);
			if(percent >= maxPercent && percent >= bestPercent){
				bestIndex = i - 1;
				bestPercent = percent;
			}
		}
		
		if(bestIndex != -1){
			s = s.substring(0, bestIndex + 1);
			q = q.substring(0, bestIndex + 1);
		}
		
		return new String[]{s, q};
	}
}

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Random;

public class UtilMethods {
	private static HashMap<Character, Character> complements = new HashMap<Character, Character>();
	private static char[] bp = {'A', 'T', 'C', 'G', 'N'};
	
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
		}
		e.printStackTrace(System.err);
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
	
	//picks a random quality value (higher)
	public static char randQualityChar(Random r){
		return (char)('0' + r.nextInt('A' - '0' + 1));
	}
	
	//compress base pairs into 3 bits each
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
	
	//decompress bits into base pairs
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
	
	//generate masks for patterns that will not change
	//this can be done once and the pattern can be used for every single text that is searched
	public static HashMap<Character, BitVector> genPatternMasks(String b, boolean indels, boolean wildcard){
		HashMap<Character, BitVector> res = new HashMap<Character, BitVector>();
		for(int i = 0; i < bp.length; i++){
			if(wildcard && bp[i] == 'N'){
				res.put(bp[i], new BitVector(b.length() + (indels ? 0 : 1)).set(0, b.length()));
			}else{
				res.put(bp[i], new BitVector(b.length() + (indels ? 0 : 1)));
			}
		}
		for(int i = 0; i < b.length(); i++){
			char c = Character.toUpperCase(b.charAt(i));
			if(wildcard && c == 'N'){
				for(int j = 0; j < bp.length; j++){
					res.get(bp[j]).set(i);
				}
			}else{
				res.get(c).set(i);
			}
		}
		return res;
	}
	
	//fuzzy search for string b in string a
	//returns a list of string ending positions and other information
	//set minOverlap to Integer.MAX_VALUE to make sure that the match only appears within the string to be searched
	//supports insertions, deletions, and substitutions, or just substitutions only
	public static ArrayList<Match> searchWithN(String a, int s, int e, String b, double edit, boolean indels, boolean bestOnly, int minOverlap, boolean wildcard, HashMap<Character, BitVector> pm){
		if(b.isEmpty())
			return new ArrayList<Match>(Arrays.asList(new Match(0, 0, 0)));
		
		minOverlap = Math.min(minOverlap, b.length());
		a = makeStr('#', b.length() - minOverlap) + a;
		e += b.length() - minOverlap;
		int min = Integer.MAX_VALUE;
		ArrayList<Match> res = new ArrayList<Match>();
		
		if(indels){ //Myer's algorithm
			BitVector vn = new BitVector(b.length());
			BitVector vp = new BitVector(b.length()).set(0, b.length());
			BitVector allSet = new BitVector(b.length()).set(0, b.length());
			int dist = b.length();
			for(int i = s; i < e; i++){
				BitVector m = a.charAt(i) == '#' ? allSet : pm.get(Character.toUpperCase(a.charAt(i)));
				BitVector d0 = new BitVector(b.length()).or(m).and(vp).add(vp).xor(vp).or(m).or(vn);
				BitVector hp = new BitVector(b.length()).or(d0).or(vp).not().or(vn);
				BitVector hn = new BitVector(b.length()).or(vp).and(d0);
				vp = new BitVector(b.length()).or(d0).orLShift(hp).not().orLShift(hn);
				vn = new BitVector(b.length()).or(d0).andLShift(hp);
				
				if(hp.get(b.length() - 1)){
					dist++;
				}else if(hn.get(b.length() - 1)){
					dist--;
				}
				
				int index = i - (b.length() - minOverlap);
				int length = Math.min(index + 1, b.length());
				if(dist <= (edit < 0.0 ? (-edit * length) : edit) && length >= minOverlap){
					if(!bestOnly || dist <= min){
						res.add(new Match(index, dist, length));
						min = dist;
					}
				}
			}
		}else{ //Bitap algorithm
			if(e - s < b.length())
				return new ArrayList<Match>();
			
			int totalEdit = (int)(edit < 0.0 ? (-edit * b.length()) : edit);
			BitVector[] r = new BitVector[totalEdit + 1];
			for(int i = 0; i <= totalEdit; i++){
				r[i] = new BitVector(b.length() + 1).set(0);
			}
			
			for(int i = s; i < e; i++){
				BitVector old = new BitVector(b.length() + 1).or(r[0]);
				boolean found = false;
				for(int j = 0; j <= totalEdit; j++){
					if(j == 0){
						if(a.charAt(i) != '#')
							r[0].and(pm.get(Character.toUpperCase(a.charAt(i))));
					}else{
						BitVector temp = new BitVector(b.length() + 1).or(r[j]);
						(a.charAt(i) == '#' ? r[j] : r[j].and(pm.get(Character.toUpperCase(a.charAt(i))))).or(old);
						old = temp;
					}
					r[j].leftShift().set(0);
					
					if(!found && r[j].get(b.length())){
						int index = i - (b.length() - minOverlap);
						int length = Math.min(index + 1, b.length());
						if(j <= (edit < 0.0 ? (-edit * length) : edit) && length >= minOverlap){
							if(!bestOnly || j <= min){
								res.add(new Match(index, j, length));
								min = j;
							}
						}
						found = true;
					}
				}
			}
		}
		
		if(bestOnly && !res.isEmpty()){ //if only find best match, then return only the matches with shortest distance
			ArrayList<Match> res2 = new ArrayList<Match>();
			for(int i = res.size() - 1; i >= 0; i--){
				if(res.get(i).edits == min){
					res2.add(res.get(i));
				}else{
					break;
				}
			}
			return res2;
		}
		return res;
	}
	
	//posterior probability based matching
	//supports one or two sequences with probability information
	//finds the best match based on the highest probability of a match that is better than the random model
	public static ArrayList<Match> searchWithProb(String a, int s, int e, String qA, String b, String qB, double prior, int minOverlap, boolean wildcard){
		if(b.isEmpty())
			return new ArrayList<Match>(Arrays.asList(new Match(0, 0, 0)));
		
		double bestProb = 0.0;
		Match bestMatch = null;
		for(int i = s + Math.max(Math.min(minOverlap, b.length()) - 1, 0); i < e; i++){
			double likelihoodRandom = 1.0;
			double likelihoodMatch = 1.0;
			for(int j = Math.max(i - b.length() + 1, 0); j <= i; j++){
				char x = Character.toUpperCase(a.charAt(j));
				char y = Character.toUpperCase(b.charAt(b.length() - 1 - (i - j)));
				double ex = toError(qA.charAt(j));
				double ey = qB == null ? 0.0 : toError(qB.charAt(b.length() - 1 - (i - j)));
				if(qB == null){
					if(!wildcard || (x != 'N' && y != 'N')){
						likelihoodMatch *= x == y ? (1 - ex) : ex;
						likelihoodRandom *= x == y ? 0.25 : (1 - 0.25);
					}
				}else{
					if(!wildcard || (x != 'N' && y != 'N')){
						likelihoodMatch *= x == y ? ((1 - ex) * (1 - ey)) : (1 - (1 - ex) * (1 - ey));
						likelihoodRandom *= x == y ? (4 * 0.25 * 0.25) : (1 - 4 * 0.25 * 0.25);
					}
				}
			}
			double total = likelihoodMatch * prior + likelihoodRandom * (1 - prior);
			double probMatch = likelihoodMatch * prior / total;
			double probRandom = likelihoodRandom * (1 - prior) / total;
			if(probMatch > probRandom && probMatch > bestProb){
				bestProb = probMatch;
				bestMatch = new Match(i, 0, Math.min(i + 1, b.length()));
			}
		}
		
		if(bestMatch == null)
			return new ArrayList<Match>();
		return new ArrayList<Match>(Arrays.asList(bestMatch));
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
	public static String[] mergeReads(String s1, String q1, String s2, String q2, double editMax, double prob, int minOverlap, boolean wildcard){
		s2 = reverseComplement(s2);
		q2 = reverse(q2);
		
		int start = -1;
		
		if(prob < 0.0){ //find match with least edit distance
			ArrayList<Match> matches = searchWithN(s2, 0, s2.length(), s1, editMax, false, false, minOverlap, wildcard, genPatternMasks(s1, false, wildcard));
			int minEdits = Integer.MAX_VALUE;
			int maxLength = 0;
			
			for(int i = 0; i < matches.size(); i++){
				if(matches.get(i).length >= maxLength && (matches.get(i).length > maxLength || matches.get(i).edits < minEdits)){
					start = s1.length() - 1 - matches.get(i).end;
					maxLength = matches.get(i).length;
					minEdits = matches.get(i).edits;
				}
			}
		}else{ //probability based matching
			ArrayList<Match> matches = searchWithProb(s2, 0, s2.length(), q2, s1, q1, prob, minOverlap, wildcard);
			if(!matches.isEmpty())
				start = s1.length() - 1 - matches.get(0).end;
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
					if((wildcard && (Character.toUpperCase(s1.charAt(i)) == 'N' || Character.toUpperCase(s2.charAt(i - start)) == 'N')) || Character.toUpperCase(s1.charAt(i)) == Character.toUpperCase(s2.charAt(i - start))){
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
	public static String[] removeAdapters(String s, String q, ArrayList<Adapter> adapters, double editMax, int minOverlap, int maxOffset, boolean indels, double prob, boolean wildcard, ArrayList<HashMap<Character, BitVector>> pm){
		int bestLength = 0;
		int bestEdit = Integer.MAX_VALUE;
		boolean bestStart = false;
		Match bestMatch = null;
		
		for(int i = 0; i < adapters.size(); i++){
			Adapter a = adapters.get(i);
			
			ArrayList<Match> matches;
			if(a.anchored){
				if(prob < 0.0)
					matches = searchWithN(a.isStart ? s : reverse(s), 0, Math.min(a.str.length() + (prob < 0.0 && indels ? (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax) : 0), s.length()), a.isStart ? a.str : reverse(a.str), editMax, indels, true, Integer.MAX_VALUE, wildcard, pm.get(i));
				else
					matches = searchWithProb(a.isStart ? s : reverse(s), 0, Math.min(a.str.length() + (prob < 0.0 && indels ? (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax) : 0), s.length()), a.isStart ? q : reverse(q), a.isStart ? a.str : reverse(a.str), null, prob, Integer.MAX_VALUE, wildcard);
			}else{
				ArrayList<Match> tempMatches = null;
				if(prob < 0.0)
					tempMatches = searchWithN(a.isStart ? s : reverse(s), 0, Math.min(maxOffset + a.str.length() + (prob < 0.0 && indels ? (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax) : 0), s.length()), a.isStart ? a.str : reverse(a.str), editMax, indels, false, minOverlap, wildcard, pm.get(i));
				else
					tempMatches = searchWithProb(a.isStart ? s : reverse(s), 0, Math.min(maxOffset + a.str.length() + (prob < 0.0 && indels ? (editMax < 0.0 ? (int)(-editMax * a.str.length()) : (int)editMax) : 0), s.length()), a.isStart ? q : reverse(q), a.isStart ? a.str : reverse(a.str), null, prob, minOverlap, wildcard);
				matches = new ArrayList<Match>();
				
				int minEdit = Integer.MAX_VALUE;
				int maxLength = 0;
				int matchIndex = -1;
				for(int j = 0; j < tempMatches.size(); j++){
					Match m = tempMatches.get(j);
					if(m.length >= maxLength && (m.length > maxLength || m.edits < minEdit)){
						matchIndex = j;
						maxLength = m.length;
						minEdit = m.edits;
					}
				}
				if(matchIndex != -1){
					matches.add(tempMatches.get(matchIndex));
				}
			}
			if(!matches.isEmpty()){
				if(matches.get(0).length >= bestLength && (matches.get(0).length > bestLength || matches.get(0).edits < bestEdit)){
					bestStart = a.isStart;
					bestMatch = matches.get(0);
					bestLength = matches.get(0).length;
					bestEdit = matches.get(0).edits;
				}
			}
		}
		
		if(bestMatch != null){
			if(bestStart){
				s = s.substring(bestMatch.end + 1);
				q = q.substring(bestMatch.end + 1);
			}else{
				s = s.substring(0, s.length() - 1 - bestMatch.end); //reverse the index to get the correct index
				q = q.substring(0, q.length() - 1 - bestMatch.end);
			}
		}
		
		return new String[]{s, q};
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
	
	public static String randSeq(Random r, int length){
		char[] res = new char[length];
		for(int i = 0; i < length; i++){
			res[i] = bp[r.nextInt(bp.length - 1)]; //no 'N' base pairs!
		}
		return new String(res);
	}
	
	public static String randQuality(Random r, int length){
		char[] res = new char[length];
		for(int i = 0; i < length; i++){
			res[i] = randQualityChar(r);
		}
		return new String(res);
	}
	
	public static String randEdit(Random r, String s, int edit, boolean indels){
		StringBuilder b = new StringBuilder(s);
		for(int i = 0; i < edit; i++){
			int idx = r.nextInt(b.length());
			int mode = indels ? r.nextInt(3) : 0;
			
			if(mode == 0){
				b.setCharAt(idx, bp[r.nextInt(bp.length - 1)]);
			}else if(mode == 1){
				b.insert(idx, bp[r.nextInt(bp.length - 1)]);
			}else{
				b.deleteCharAt(idx);
			}
		}
		
		return b.toString();
	}
}

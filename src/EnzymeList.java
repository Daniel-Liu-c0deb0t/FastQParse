import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class EnzymeList {
	public static final HashMap<String, ArrayList<String>> enzymes = new HashMap<String, ArrayList<String>>();
	
	public static void init(){
		enzymes.put("APEKI", new ArrayList<String>(Arrays.asList("CAGC", "CTGC")));
		enzymes.put("PSTI", new ArrayList<String>(Arrays.asList("TGCAG")));
		enzymes.put("ECOT22I", new ArrayList<String>(Arrays.asList("TGCAT")));
		enzymes.put("PASI", new ArrayList<String>(Arrays.asList("CAGGG", "CTGGG")));
		enzymes.put("HPAII", new ArrayList<String>(Arrays.asList("CGG")));
		enzymes.put("MSPI", new ArrayList<String>(Arrays.asList("CGG")));
		enzymes.put("PSTI-ECOT22I", new ArrayList<String>(Arrays.asList("TGCAG", "TGCAT")));
		enzymes.put("PSTI-MSPI", new ArrayList<String>(Arrays.asList("TGCAG")));
		enzymes.put("PSTI-TAQI", new ArrayList<String>(Arrays.asList("TGCAG")));
		enzymes.put("SBFI-MSPI", new ArrayList<String>(Arrays.asList("TGCAGG")));
		enzymes.put("ASISI-MSPI", new ArrayList<String>(Arrays.asList("ATCGC")));
		enzymes.put("BSSHII-MSPI", new ArrayList<String>(Arrays.asList("CGCGC")));
		enzymes.put("FSEI-MSPI", new ArrayList<String>(Arrays.asList("CCGGCC")));
		enzymes.put("SALI-MSPI", new ArrayList<String>(Arrays.asList("TCGAC")));
		enzymes.put("APOI", new ArrayList<String>(Arrays.asList("AATTC", "AATTT")));
		enzymes.put("BAMHI", new ArrayList<String>(Arrays.asList("GATCC")));
		enzymes.put("MSEI", new ArrayList<String>(Arrays.asList("TAA")));
		enzymes.put("SAU3AI", new ArrayList<String>(Arrays.asList("GATC")));
		enzymes.put("RBSTA", new ArrayList<String>(Arrays.asList("TA")));
		enzymes.put("RBSCG", new ArrayList<String>(Arrays.asList("CG")));
		enzymes.put("NSPI", new ArrayList<String>(Arrays.asList("CATGT", "CATGC")));
		enzymes.put("AVAII", new ArrayList<String>(Arrays.asList("GACC", "GTCC")));
		enzymes.put("NA", new ArrayList<String>(Arrays.asList("")));
	}
}

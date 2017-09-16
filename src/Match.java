
public class Match{
	public int end, edits, correctLength;
	
	public Match(int end, int edits, int correctLength){
		this.end = end;
		this.edits = edits;
		this.correctLength = correctLength;
	}
	
	@Override
	public String toString(){
		return "End Index: " + end + "; Edits: " + edits + "; Correct Length: " + correctLength;
	}
}

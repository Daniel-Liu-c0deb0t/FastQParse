
public class Match{
	public int start, edits, correctLength;
	
	public Match(int start, int edits, int correctLength){
		this.start = start;
		this.edits = edits;
		this.correctLength = correctLength;
	}
	
	@Override
	public String toString(){
		return "Start Index: " + start + "; Edits: " + edits + "; Correct Length: " + correctLength;
	}
}

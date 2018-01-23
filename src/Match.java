
public class Match{
	public int end, edits, length;
	
	//note that the length is not the exact length of the match!
	public Match(int end, int edits, int length){
		this.end = end;
		this.edits = edits;
		this.length = length;
	}
	
	@Override
	public String toString(){
		return "End Index: " + end + "; Edits: " + edits + "; Length: " + length;
	}
}

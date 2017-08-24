
public class Adapter {
	public String str;
	public boolean isStart;
	public boolean anchored;
	
	public Adapter(String str, boolean isStart, boolean anchored){
		this.str = str;
		this.isStart = isStart;
		this.anchored = anchored;
	}
	
	@Override
	public String toString(){
		return (isStart ? "5' " : "3' ") + (anchored ? "anchored " : "") + str;
	}
}

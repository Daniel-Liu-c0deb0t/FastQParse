import java.util.ArrayList;

public class Strings{
	public ArrayList<String> str = new ArrayList<String>();
	
	public void add(String s){
		str.add(s);
	}
	
	public String get(int i){
		return str.get(i);
	}
	
	public int size(){
		return str.size();
	}
}

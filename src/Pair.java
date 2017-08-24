
public class Pair<T>{
	public T a, b;
	
	public Pair(T a, T b){
		this.a = a;
		this.b = b;
	}
	
	@Override
	public String toString(){
		return a.toString() + " " + b.toString();
	}
}

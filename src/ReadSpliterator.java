import java.io.BufferedReader;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Consumer;

public class ReadSpliterator<T> implements Spliterator<T>{
	private int batch;
	private BufferedReader rF, rR, rIF, rIR;
	
	//ReadSpliterator's type is needs to be Read, or else bad stuff will happen
	public ReadSpliterator(int batch, BufferedReader rF, BufferedReader rR, BufferedReader rIF, BufferedReader rIR){
		this.batch = batch;
		this.rF = rF;
		this.rR = rR;
		this.rIF = rIF;
		this.rIR = rIR;
	}
	
	@Override
	public int characteristics(){
		return NONNULL | CONCURRENT | IMMUTABLE;
	}
	
	@Override
	public long estimateSize(){
		return Long.MAX_VALUE;
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public boolean tryAdvance(Consumer<? super T> c){
		String[] nextF = new String[4];
		String[] nextR = null;
		String[] nextIF = null;
		String[] nextIR = null;
		
		try{
			for(int i = 0; i < 4; i++)
				nextF[i] = rF.readLine();
			if(rR != null){
				nextR = new String[4];
				for(int i = 0; i < 4; i++)
					nextR[i] = rR.readLine();
			}
			if(rIF != null){
				nextIF = new String[4];
				for(int i = 0; i < 4; i++)
					nextIF[i] = rIF.readLine();
			}
			if(rIR != null){
				nextIR = new String[4];
				for(int i = 0; i < 4; i++)
					nextIR[i] = rIR.readLine();
			}
			if(nextF[3] == null || (rR != null && nextR[3] == null) || (rIF != null && nextIF[3] == null) || (rIR != null && nextIR[3] == null)){
				return false;
			}
		}catch(Exception e){
			UtilMethods.defaultExceptionHandler(null, e);
		}
		
		c.accept((T)new Read(nextF, nextR, nextIF, nextIR));
		return true;
	}
	
	@Override
	public Spliterator<T> trySplit(){
		TempHolder<T> h = new TempHolder<T>();
		if(!tryAdvance(h)){
			return null;
		}
		Read[] o = new Read[batch];
		int i = 0;
		do{
			o[i] = h.get();
		}while(++i < batch && tryAdvance(h));
		return Spliterators.spliterator(o, 0, i, characteristics() | SIZED);
	}
	
	@SuppressWarnings("hiding")
	private class TempHolder<T> implements Consumer<T>{
		private Read val;
		
		@Override
		public void accept(T val){
			this.val = (Read)val;
		}
		
		public Read get(){
			return val;
		}
	}
}

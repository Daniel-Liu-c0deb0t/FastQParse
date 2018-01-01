import java.util.BitSet;

public class CompressedRead{
	public double error;
	
	public String descriptionFRead;
	public BitSet sequenceFRead;
	public byte[] qualityFRead;
	
	public String descriptionRRead;
	public BitSet sequenceRRead;
	public byte[] qualityRRead;
	
	public String descriptionFIndex;
	public BitSet sequenceFIndex;
	public byte[] qualityFIndex;
	
	public String descriptionRIndex;
	public BitSet sequenceRIndex;
	public byte[] qualityRIndex;
	
	public CompressedRead(double error, String descriptionFRead, String sequenceFRead, String qualityFRead,
			String descriptionFIndex, String sequenceFIndex, String qualityFIndex){
		this.error = error;
		
		this.descriptionFRead = descriptionFRead;
		this.sequenceFRead = UtilMethods.toBit(sequenceFRead);
		this.qualityFRead = UtilMethods.qScoreToByteArray(qualityFRead);
		
		this.descriptionFIndex = descriptionFIndex;
		this.sequenceFIndex = UtilMethods.toBit(sequenceFIndex);
		this.qualityFIndex = UtilMethods.qScoreToByteArray(qualityFIndex);
	}
	
	public CompressedRead(double error, String descriptionFRead, String sequenceFRead, String qualityFRead,
			String descriptionFIndex, String sequenceFIndex, String qualityFIndex,
			String descriptionRRead, String sequenceRRead, String qualityRRead,
			String descriptionRIndex, String sequenceRIndex, String qualityRIndex){
		this.error = error;
		
		this.descriptionFRead = descriptionFRead;
		this.sequenceFRead = UtilMethods.toBit(sequenceFRead);
		this.qualityFRead = UtilMethods.qScoreToByteArray(qualityFRead);
		
		this.descriptionFIndex = descriptionFIndex;
		this.sequenceFIndex = UtilMethods.toBit(sequenceFIndex);
		this.qualityFIndex = UtilMethods.qScoreToByteArray(qualityFIndex);
		
		this.descriptionRRead = descriptionRRead;
		this.sequenceRRead = UtilMethods.toBit(sequenceRRead);
		this.qualityRRead = UtilMethods.qScoreToByteArray(qualityRRead);
		
		this.descriptionRIndex = descriptionRIndex;
		this.sequenceRIndex = UtilMethods.toBit(sequenceRIndex);
		this.qualityRIndex = UtilMethods.qScoreToByteArray(qualityRIndex);
	}
}

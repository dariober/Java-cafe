package samTextViewer;

public class TextWindow {
	
	// Genomic coords spanend by window
	private String chrom;
	private int from;
	private int to;
	
	/*                  Constructors                 */
	public TextWindow() {}
	
	public TextWindow(int from, int to) {
		if(from > to){
			System.err.println("Invalid window range: from > to"); 
			System.exit(1);
		}
		this.from= from;
		this.to= to;
	}
	
	/* Methods */
	
	public String toString(){
		String s= chrom + ":" + from + "-" + to;
		return s;
	}
	
	/*                  Setters & Getters             */
	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public int getFrom() {
		return from;
	}

	public void setFrom(int from) {
		this.from = from;
	}

	public int getTo() {
		return to;
	}

	public void setTo(int to) {
		this.to = to;
	}
	
}

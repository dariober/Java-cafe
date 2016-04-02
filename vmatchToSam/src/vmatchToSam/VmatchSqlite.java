package vmatchToSam;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class VmatchSqlite {

	private File sqlitedb= null; 
	private Connection conn= null;
	private String sqlInsert= "INSERT INTO " + VmatchRecord.SQLITE_VMATCH_TABLE + " ("
			+ "referenceAlignmentLength, "
			+ "referenceNameOrIndex, "
			+ "referenceAlignmentStart, "
			+ "strand, "
			+ "queryAlignmentLength, "
			+ "queryNameOrIndex, "
			+ "queryAlignmentStart, "
			+ "editDistance, "
			+ "evalue, "
			+ "alignmentScore, "
			+ "pctIdentity, "
			+ "alignedQuerySequence, "
			+ "alignedReferenceSequence"
			+ ") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";

	/* G E T T E R S   AND   S E T T E R S */
	protected File getSqlitedb() {
		return sqlitedb;
	}

	protected String getSqlInsert() {
		return sqlInsert;
	}

	
	protected Connection getConn() {
		return conn;
	}
	
	/* C O N S T R U C T O R */
	public VmatchSqlite() throws ClassNotFoundException, IOException, SQLException{
		this.createSQLiteVmatchDb();
		this.conn= DriverManager.getConnection("jdbc:sqlite:" + this.sqlitedb.getName());
	}
	
	private void createSQLiteVmatchDb() throws IOException, ClassNotFoundException, SQLException {
		
	    // Get a tmp file name: Delete tmp file, keep the name.
	    File tmpfile= File.createTempFile("vmatchToSam.", ".tmp.db");
	    String sqlitedbName= tmpfile.getName();
	    tmpfile.delete();
	    
	    File sqlitedb= new File(sqlitedbName);
	    sqlitedb.deleteOnExit();

	    Class.forName("org.sqlite.JDBC");
	    Connection conn = DriverManager.getConnection("jdbc:sqlite:" + sqlitedb.getName());
        Statement stmt = conn.createStatement();
        String sql = "CREATE TABLE " +
                   VmatchRecord.SQLITE_VMATCH_TABLE + " (" +
                   "referenceAlignmentLength int, " +
                   "referenceNameOrIndex text, " +
                   "referenceAlignmentStart int, " + 
                   "strand char, " +
                   "queryAlignmentLength int, " +
                   "queryNameOrIndex text, " +
                   "queryAlignmentStart int, " +
                   "editDistance int, " +
                   "evalue float, " +
                   "alignmentScore int, " +
                   "pctIdentity float, " +
                   "alignedQuerySequence text, " +
                   "alignedReferenceSequence text" +
                   ")"; 
        stmt.executeUpdate(sql);
        stmt.close();   
        conn.close();
        this.sqlitedb= sqlitedb;
	}
	
	public void add(VmatchRecord vmatchRecord, PreparedStatement stmt) throws SQLException {
		// It's important to get the prepared statement from outside so it is not recompiled
		// for each execution of .add().
		stmt.setInt(1,     vmatchRecord.getReferenceAlignmentLength());
		stmt.setString(2,  vmatchRecord.getReferenceNameOrIndex());
		stmt.setInt(3,     vmatchRecord.getQueryAlignmentStart());
		stmt.setString(4,  vmatchRecord.getStrand());
		stmt.setInt(5,     vmatchRecord.getQueryAlignmentLength());
		stmt.setString(6,  vmatchRecord.getQueryNameOrIndex());
		stmt.setInt(7,     vmatchRecord.getQueryAlignmentStart());
		stmt.setInt(8,     vmatchRecord.getEditDistance());
		stmt.setFloat(9,   vmatchRecord.getEvalue());
		stmt.setInt(10,    vmatchRecord.getAlignmentScore());
		stmt.setFloat(11,  vmatchRecord.getPctIdentity());
		stmt.setString(12, vmatchRecord.getAlignedQuerySequence());
		stmt.setString(13, vmatchRecord.getAlignedReferenceSequence());
	}
}

package vmatchToSam;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import org.junit.Test;

import exceptions.InvalidVmatchRecordException;

public class VmatchSqliteTest {
	
	@Test
	public void canCreateDatabase() throws IOException, SQLException, ClassNotFoundException { 

		VmatchSqlite vsql= new VmatchSqlite(); 
		assertTrue(vsql.getSqlitedb().exists());
		
	}
	
	@Test
	public void canAddVmatchAlignmentToDb() throws IOException, InvalidVmatchRecordException, SQLException, ClassNotFoundException, InstantiationException, IllegalAccessException{

		BufferedReader br= new BufferedReader(new FileReader(new File("test_data/aln_1.vmatch.txt")));
		VmatchRecord vRec= new VmatchRecord(br);
		
		VmatchSqlite vsql= new VmatchSqlite();
		PreparedStatement stmtInsert= vsql.getConn().prepareStatement(vsql.getSqlInsert());
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
	
	}
	
	@Test
	public void canRetrieveResults() throws IOException, InvalidVmatchRecordException, ClassNotFoundException, SQLException{

		BufferedReader br= new BufferedReader(new FileReader(new File("test_data/aln_1.vmatch.txt")));
		VmatchRecord vRec= new VmatchRecord(br);
		
		VmatchSqlite vsql= new VmatchSqlite();
		PreparedStatement stmtInsert= vsql.getConn().prepareStatement(vsql.getSqlInsert());
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		vsql.add(vRec, stmtInsert); stmtInsert.execute();
		
		Statement stmt= vsql.getConn().createStatement();
		ResultSet rs= stmt.executeQuery(
				"select * from " + VmatchRecord.SQLITE_VMATCH_TABLE + 
				" order by queryNameOrIndex, alignmentScore");
	
		while(rs.next()){
			VmatchRecord vRecSql= new VmatchRecord(rs);
			assertEquals(vRec.getQueryNameOrIndex(), vRecSql.getQueryNameOrIndex());
		}
	}
}

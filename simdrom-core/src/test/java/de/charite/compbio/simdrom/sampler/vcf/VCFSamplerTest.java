/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.base.Optional;

import de.charite.compbio.simdrom.interval.SAMFileHeaderBuilder;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFSamplerTest {

	private static File file;
	private static VCFFileReader reader;
	private VCFSampler samplerAF;
	private VCFSampler samplerHasNext;
	private VCFSampler samplerNoNext;
	private VCFSampler samplerProb;
	private VCFSampler samplerAmount;
	private VCFSampler samplerAC_AN;
	private VCFSampler samplerProb_interval;

	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		file = new File(VCFSamplerTest.class.getResource("/ExAC.r0.3.sites.vep.head300.vcf.gz").toURI().getPath());
		reader = new VCFFileReader(file);
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		samplerHasNext = new VCFSampler.Builder().probability(1.0).seed(42).vcfReader(reader).build();
		samplerNoNext = new VCFSampler.Builder().probability(0.0).seed(42).vcfReader(reader).build();
		samplerAF = new VCFSampler.Builder().afIdentifier("AF").seed(42).vcfReader(reader).build();
		samplerProb = new VCFSampler.Builder().probability(0.5).seed(42).vcfReader(reader).build();
		samplerAmount = new VCFSampler.Builder().variantsAmount(10).seed(42).vcfReader(reader).build();
		samplerAC_AN = new VCFSampler.Builder().acIdentifier("AC").anIdentifier("AN").seed(42).vcfReader(reader)
				.build();
		IntervalList list = new IntervalList(SAMFileHeaderBuilder.build());
		list.add(new Interval("1", 69649, 69655));
		samplerProb_interval = new VCFSampler.Builder().probability(1.0).seed(42).intervals(list).vcfReader(reader).build();
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#hasNext()}.
	 */
	@Test
	public void testHasNext() {
		assertTrue(samplerHasNext.hasNext());
		assertFalse(samplerNoNext.hasNext());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#next()}.
	 */
	@Test
	public void testNext() {
		int count = 0;
		while (samplerProb.hasNext()) {
			samplerProb.next();
			count++;
		}
		assertEquals(53, count);
		count = 0;
		while (samplerProb_interval.hasNext()) {
			samplerProb_interval.next();
			count++;
		}
		assertEquals(2, count);
		count = 0;
		while (samplerAmount.hasNext()) {
			samplerAmount.next();
			count++;
		}
		assertEquals(10, count);
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getSample()}.
	 */
	@Test
	public void testGetSample() {
		assertFalse(samplerHasNext.getSample().isPresent());
		// TODO cehck real sample
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getProbability()}.
	 */
	@Test
	public void testGetProbability() {
		assertEquals(0.5, samplerProb.getProbability(),0.0);
		assertEquals(1.0, samplerHasNext.getProbability(),0.0);
		assertEquals(0.0, samplerNoNext.getProbability(),0.0);
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getVCFHeader()}.
	 */
	@Test
	public void testGetVCFHeader() {
		fail("Not yet implemented");
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getSampleNames()}.
	 */
	@Test
	public void testGetSampleNames() {
		assertEquals(1, samplerHasNext.getSampleNames().size());
		assertTrue(samplerHasNext.getSampleNames().contains("Sampled"));
		
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getVariantsAmount()}.
	 */
	@Test
	public void testGetVariantsAmount() {
		assertEquals(0, samplerHasNext.getVariantsAmount());
		assertEquals(10, samplerAmount.getVariantsAmount());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getAFIdentifier()}.
	 */
	@Test
	public void testGetAFIdentifier() {
		assertEquals("AF", samplerAF.getAFIdentifier());
		assertNull(samplerAC_AN.getAFIdentifier());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getFilters()}.
	 */
	@Test
	public void testGetFilters() {
		fail("Not yet implemented");
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getACIdentifier()}.
	 */
	@Test
	public void testGetACIdentifier() {
		assertEquals("AC", samplerAC_AN.getACIdentifier());
		assertNull(samplerAF.getACIdentifier());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getANIdentifier()}.
	 */
	@Test
	public void testGetANIdentifier() {
		assertEquals("AN", samplerAC_AN.getANIdentifier());
		assertNull(samplerAF.getANIdentifier());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getIntervals()}.
	 */
	@Test
	public void testGetIntervals() {
		assertEquals(1, samplerProb_interval.getIntervals().size());
		assertEquals(0, samplerAF.getIntervals().size());
	}

}

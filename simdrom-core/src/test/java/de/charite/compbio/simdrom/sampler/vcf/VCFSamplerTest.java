/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.base.Charsets;
import com.google.common.collect.ImmutableSet;
import com.google.common.io.Files;

import de.charite.compbio.simdrom.exception.InfoIDNotFoundException;
import de.charite.compbio.simdrom.filter.GreaterOrEqualInfoFieldFilter;
import de.charite.compbio.simdrom.interval.SAMFileHeaderBuilder;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFSamplerTest {

	private static File exacFile;
	private static File kgFile;
	private static File headerSampler;
	private static VCFFileReader exacReader;
	private static VCFFileReader kgReader;
	private VCFSampler samplerAF;
	private VCFSampler samplerAF1KG;
	private VCFSampler samplerHasNext;
	private VCFSampler samplerNoNext;
	private VCFSampler samplerProb;
	private VCFSampler samplerFilter;
	private VCFSampler samplerSample;
	private VCFSampler samplerAmount;
	private VCFSampler samplerAC_AN;
	private VCFSampler samplerProb_interval;

	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		exacFile = new File(VCFSamplerTest.class.getResource("/ExAC.r0.3.sites.vep.head300.vcf.gz").toURI().getPath());
		kgFile = new File(VCFSamplerTest.class
				.getResource("/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.head300.vcf.gz")
				.toURI().getPath());
		exacReader = new VCFFileReader(exacFile);
		kgReader = new VCFFileReader(kgFile);

		headerSampler = new File(VCFSamplerTest.class.getResource("/Exac.header.txt").toURI().getPath());
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {

		samplerHasNext = new VCFSampler.Builder().probability(1.0).seed(42).vcfReader(exacReader).build();
		samplerNoNext = new VCFSampler.Builder().probability(0.0).seed(42).vcfReader(exacReader).build();
		samplerAF1KG = new VCFSampler.Builder().afIdentifier("AF").seed(42).vcfReader(kgReader).build();
		samplerAF = new VCFSampler.Builder().afIdentifier("AF").seed(42).vcfReader(exacReader).build();
		samplerSample = new VCFSampler.Builder().sample("HG00097").seed(42).vcfReader(kgReader).build();
		samplerFilter = new VCFSampler.Builder().sample("HG00097")
				.filters(ImmutableSet.of(new GreaterOrEqualInfoFieldFilter("AF", 0.5))).seed(42).vcfReader(kgReader)
				.build();
		samplerProb = new VCFSampler.Builder().probability(0.5).seed(42).vcfReader(exacReader).build();
		samplerAmount = new VCFSampler.Builder().variantsAmount(10).seed(42).vcfReader(exacReader).build();
		samplerAC_AN = new VCFSampler.Builder().acIdentifier("AC").anIdentifier("AN").seed(42).vcfReader(exacReader)
				.build();
		IntervalList list = new IntervalList(SAMFileHeaderBuilder.build());
		list.add(new Interval("1", 69649, 69655));
		samplerProb_interval = new VCFSampler.Builder().probability(1.0).seed(42).intervals(list).vcfReader(exacReader)
				.build();
	}

	/**
	 * Test {@link VCFSampler.Builder}
	 */
	@Test(expected = InfoIDNotFoundException.class)
	public void testBuilder() {
		new VCFSampler.Builder().afIdentifier("AlleleF").seed(42).vcfReader(exacReader).build();
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
		// assertEquals(10, count);

		count = 0;
		while (samplerAF1KG.hasNext()) {
			samplerAF1KG.next();
			count++;
		}
		assertEquals(3, count);

		count = 0;
		while (samplerAF.hasNext()) {
			samplerAF.next();
			count++;
		}
		assertEquals(1, count);

		count = 0;
		while (samplerAC_AN.hasNext()) {
			samplerAC_AN.next();
			count++;
		}
		assertEquals(1, count);

		count = 0;
		while (samplerSample.hasNext()) {
			samplerSample.next();
			count++;
		}
		assertEquals(9, count);

		count = 0;
		while (samplerFilter.hasNext()) {
			samplerFilter.next();
			count++;
		}
		assertEquals(1, count);
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getSample()}.
	 */
	@Test
	public void testGetSample() {
		assertFalse(samplerHasNext.getSample().isPresent());
		assertTrue(samplerSample.getSample().isPresent());
		assertEquals("HG00097", samplerSample.getSample().get());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getProbability()}.
	 */
	@Test
	public void testGetProbability() {
		assertEquals(0.5, samplerProb.getProbability(), 0.0);
		assertEquals(1.0, samplerHasNext.getProbability(), 0.0);
		assertEquals(0.0, samplerNoNext.getProbability(), 0.0);
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getVCFHeader()}.
	 * 
	 * @throws IOException
	 */
	@Test
	public void testGetVCFHeader() throws IOException {
		assertEquals(Files.toString(headerSampler, Charsets.UTF_8), samplerHasNext.getVCFHeader().toString());
	}

	/**
	 * Test method for {@link de.charite.compbio.simdrom.sampler.vcf.VCFSampler#getSampleNames()}.
	 */
	@Test
	public void testGetSampleNames() {
		assertEquals(1, samplerHasNext.getSampleNames().size());
		assertTrue(samplerHasNext.getSampleNames().contains("Sampled"));

		assertEquals(1, samplerSample.getSampleNames().size());
		assertTrue(samplerSample.getSampleNames().contains("HG00097"));

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

		assertEquals(1, samplerFilter.getFilters().size());
		assertTrue(samplerAC_AN.getFilters().isEmpty());

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

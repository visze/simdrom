/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.Random;

import htsjdk.variant.vcf.VCFFileReader;

/**
 * This parser selects one random sample out of a VCF file. If the sample is set it only checks if it is present.
 * 
 * @author Max Schubach {@literal <max.schubach@charite.de>}
 *
 */
public class VCFSampleSelector {

	/**
	 * Seed for randomly picking up a sample. Important for debugging.
	 */
	private long seed;
	/**
	 * The selected sample. {@link Optional#isPresent()} will be <code>false</code> if a random sample is not selected.
	 */
	private Optional<String> sample;
	private List<String> samples;

	private VCFSampleSelector(File file) {
		samples = getVCFSamples(file);
	}

	/**
	 * Constructor if you want to select a random sample in the VCF file.
	 * 
	 * @param file
	 *            VCF file
	 * @param seed
	 *            Seed to pick a random sample
	 */
	public VCFSampleSelector(File file, long seed) {
		this(file);
		this.seed = seed;
	}

	/**
	 * 
	 * Constructor if you want to select a specific sample in the VCF file. It only checks if it is present.
	 * 
	 * @param file
	 *            VCF file
	 * @param sample
	 *            The sample that should be present in the header
	 */
	public VCFSampleSelector(File file, String sample) {
		this(file);
		this.sample = Optional.of(sample);
		if (!samples.contains(this.sample.get()))
			throw new RuntimeException(
					"Sample " + this.sample.get() + " is nor present in the VCF file " + file.getAbsolutePath());
	}

	/**
	 * Getter for the random/checked sample.
	 * 
	 * @return The sample that should be used.
	 */
	public String getSample() {
		if (!sample.isPresent())
			sample = selectSample();
		return sample.get();
	}

	private Optional<String> selectSample() {
		int num = new Random(seed).nextInt(samples.size());
		return Optional.of(samples.get(num));
	}

	/**
	 * using a vcf file to get the samples listed in the header. it will close the parser afterwards
	 * 
	 * @param file
	 *            The VCF file. Index not needed
	 * @return All genotyped samples that are in the header.
	 */
	private List<String> getVCFSamples(File file) {
		VCFFileReader parser = new VCFFileReader(file, false);
		List<String> samples = parser.getFileHeader().getGenotypeSamples();
		parser.close();
		return samples;
	}

}

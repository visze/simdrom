package de.charite.compbio.simdrom.cli.exception;


public class WrongIntervalFormatException extends Exception {
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 5514299661602182053L;

	public WrongIntervalFormatException(String input) {
		super(buildMessage(input));

		
	}

	private static String buildMessage(String input) {
		StringBuilder message = new StringBuilder("Could not parse interval ");
		message.append(input);
		message.append("! Format should be like chr2:121312-121312.");
		return message.toString();
	}

}
